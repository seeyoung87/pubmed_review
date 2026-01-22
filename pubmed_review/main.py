import email
import imaplib
import json
import os
import re
import ssl
from dataclasses import dataclass
from datetime import datetime, timedelta
from email.header import decode_header
from typing import Iterable
import xml.etree.ElementTree as ET

import requests
import yaml
from bs4 import BeautifulSoup
from google.oauth2 import service_account
from googleapiclient.discovery import build
from openai import OpenAI

PMID_REGEX = re.compile(r"PMID:\s*(\d+)", re.IGNORECASE)
PMID_URL_REGEX = re.compile(r"pubmed/(\d+)", re.IGNORECASE)
SEARCH_NAME_REGEX = re.compile(r"What's new for '(.+?)' in PubMed", re.IGNORECASE)
DEFAULT_IMAP_HOST = "imap.gmail.com"
DEFAULT_IMAP_PORT = 993
DEFAULT_PUBMED_FROM = "pubmed@ncbi.nlm.nih.gov"
DEFAULT_PUBMED_SUBJECT_TERMS = ("What's new for", "in PubMed")
DEFAULT_SHEET_NAME = "PubMed"


@dataclass
class Article:
    pmid: str
    title: str
    journal: str
    pub_date: str
    doi: str
    authors: str
    abstract: str
    url: str


@dataclass
class ReviewResult:
    article: Article
    high_if: bool
    is_novel: bool
    novelty_reason: str
    summary: str
    strengths: str


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)


def decode_mime_words(text: str) -> str:
    decoded = []
    for part, encoding in decode_header(text):
        if isinstance(part, bytes):
            decoded.append(part.decode(encoding or "utf-8", errors="ignore"))
        else:
            decoded.append(part)
    return "".join(decoded)


def extract_text_from_email(message: email.message.Message) -> str:
    body_parts = []
    for part in message.walk():
        content_type = part.get_content_type()
        if content_type == "text/plain" and not part.get("Content-Disposition"):
            charset = part.get_content_charset() or "utf-8"
            body_parts.append(part.get_payload(decode=True).decode(charset, errors="ignore"))
        if content_type == "text/html" and not part.get("Content-Disposition"):
            charset = part.get_content_charset() or "utf-8"
            html = part.get_payload(decode=True).decode(charset, errors="ignore")
            soup = BeautifulSoup(html, "html.parser")
            body_parts.append(soup.get_text("\n"))
    return "\n".join(body_parts)


def find_pmids(text: str) -> list[str]:
    pmids = set(PMID_REGEX.findall(text))
    pmids.update(PMID_URL_REGEX.findall(text))
    return sorted(pmids)


def extract_alert_section(text: str) -> str:
    marker = "PubMed Results"
    if marker not in text:
        return text
    _, tail = text.split(marker, 1)
    return f"{marker}{tail}"


def extract_search_name(subject: str) -> str:
    match = SEARCH_NAME_REGEX.search(subject)
    if match:
        return match.group(1).strip()
    return ""


def fetch_pubmed_summary(pmids: Iterable[str]) -> dict:
    if not pmids:
        return {}
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        params={"db": "pubmed", "id": ",".join(pmids), "retmode": "json"},
        timeout=30,
    )
    response.raise_for_status()
    return response.json().get("result", {})


def fetch_pubmed_abstracts(pmids: Iterable[str]) -> dict[str, str]:
    pmid_list = list(pmids)
    if not pmid_list:
        return {}
    response = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={"db": "pubmed", "id": ",".join(pmid_list), "retmode": "xml"},
        timeout=30,
    )
    response.raise_for_status()
    root = ET.fromstring(response.text)
    abstracts_by_pmid: dict[str, str] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid_element = article.find(".//MedlineCitation/PMID")
        if pmid_element is None or not pmid_element.text:
            continue
        abstracts = []
        for abstract in article.findall(".//AbstractText"):
            if abstract.text:
                abstracts.append(abstract.text.strip())
        abstracts_by_pmid[pmid_element.text.strip()] = "\n".join(abstracts)
    return abstracts_by_pmid


def build_article(pmid: str, summary: dict, abstracts: dict[str, str]) -> Article:
    item = summary.get(pmid, {})
    title = item.get("title", "").strip()
    journal = item.get("fulljournalname", "").strip()
    pub_date = item.get("pubdate", "").strip()
    doi = ""
    for article_id in item.get("articleids", []) or []:
        if article_id.get("idtype") == "doi":
            doi = article_id.get("value", "")
            break
    authors = "; ".join(author.get("name", "") for author in item.get("authors", []) or [])
    abstract = abstracts.get(pmid, "")
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    return Article(
        pmid=pmid,
        title=title,
        journal=journal,
        pub_date=pub_date,
        doi=doi,
        authors=authors,
        abstract=abstract,
        url=url,
    )


def is_high_if(journal: str, high_if_list: list[str]) -> bool:
    normalized = journal.lower()
    return any(name.lower() in normalized for name in high_if_list)


def openai_client() -> OpenAI:
    return OpenAI(api_key=os.environ["OPENAI_API_KEY"])


def normalize_env_value(value: str) -> str:
    return value.replace("\u00a0", " ").strip()


def parse_json_response(content: str, context: str) -> dict:
    try:
        return json.loads(content)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Failed to parse JSON for {context}") from exc


def llm_novelty(client: OpenAI, config: dict, article: Article) -> tuple[bool, str]:
    prompt = config["llm"]["novelty_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    response = client.chat.completions.create(
        model=config["llm"]["model"],
        temperature=config["llm"].get("temperature", 0.2),
        messages=[{"role": "user", "content": prompt}],
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "novelty_assessment",
                "schema": {
                    "type": "object",
                    "properties": {
                        "is_novel": {"type": "boolean"},
                        "reason": {"type": "string"},
                    },
                    "required": ["is_novel", "reason"],
                    "additionalProperties": False,
                },
            },
        },
        max_tokens=400,
    )
    content = response.choices[0].message.content
    data = parse_json_response(content, "novelty response")
    return bool(data.get("is_novel")), data.get("reason", "")


def llm_summary(client: OpenAI, config: dict, article: Article) -> tuple[str, str]:
    prompt = config["llm"]["summary_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    response = client.chat.completions.create(
        model=config["llm"]["model"],
        temperature=config["llm"].get("temperature", 0.2),
        messages=[{"role": "user", "content": prompt}],
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "summary_response",
                "schema": {
                    "type": "object",
                    "properties": {
                        "summary": {"type": "string"},
                        "strengths": {"type": "string"},
                    },
                    "required": ["summary", "strengths"],
                    "additionalProperties": False,
                },
            },
        },
        max_tokens=500,
    )
    content = response.choices[0].message.content
    data = parse_json_response(content, "summary response")
    return data.get("summary", ""), data.get("strengths", "")


def google_sheets_service() -> object:
    json_data = os.environ.get("GOOGLE_SERVICE_ACCOUNT_JSON")
    file_path = os.environ.get("GOOGLE_SERVICE_ACCOUNT_FILE")
    if json_data:
        info = json.loads(json_data)
    elif file_path:
        with open(file_path, "r", encoding="utf-8") as file:
            info = json.load(file)
    else:
        raise RuntimeError("Missing GOOGLE_SERVICE_ACCOUNT_JSON or GOOGLE_SERVICE_ACCOUNT_FILE")
    credentials = service_account.Credentials.from_service_account_info(
        info, scopes=["https://www.googleapis.com/auth/spreadsheets"]
    )
    return build("sheets", "v4", credentials=credentials)


def resolve_sheet_name(config: dict, search_name: str) -> str:
    if search_name:
        return search_name
    sheet_config = config.get("sheets", {})
    return sheet_config.get("sheet_name") or DEFAULT_SHEET_NAME


def append_rows(config: dict, sheet_name: str, rows: list[list[str]]) -> None:
    if not rows:
        return
    service = google_sheets_service()
    sheet_id = os.environ.get("SPREADSHEET_ID") or config["sheets"]["spreadsheet_id"]
    body = {"values": rows}
    service.spreadsheets().values().append(
        spreadsheetId=sheet_id,
        range=f"{sheet_name}!A1",
        valueInputOption="RAW",
        insertDataOption="INSERT_ROWS",
        body=body,
    ).execute()


def build_rows(results: list[ReviewResult]) -> list[list[str]]:
    rows = []
    for result in results:
        article = result.article
        selection = []
        if result.high_if:
            selection.append("High IF")
        if result.is_novel:
            selection.append("Novelty")
        rows.append(
            [
                datetime.utcnow().strftime("%Y-%m-%d"),
                article.pmid,
                article.title,
                article.journal,
                article.pub_date,
                article.doi,
                ", ".join(selection),
                result.novelty_reason,
                result.summary,
                result.strengths,
            ]
        )
    return rows


def connect_imap(config: dict) -> imaplib.IMAP4_SSL:
    context = ssl.create_default_context()
    email_config = config.get("email", {})
    client = imaplib.IMAP4_SSL(
        email_config.get("imap_host", DEFAULT_IMAP_HOST),
        email_config.get("imap_port", DEFAULT_IMAP_PORT),
        ssl_context=context,
    )
    email_user = normalize_env_value(os.environ["EMAIL_USER"])
    email_pass = normalize_env_value(os.environ["EMAIL_PASS"])
    client.login(email_user, email_pass)
    return client


def search_emails(client: imaplib.IMAP4_SSL, config: dict) -> list[bytes]:
    email_config = config.get("email", {})
    client.select(email_config.get("mailbox", "INBOX"))
    since_date = (datetime.utcnow() - timedelta(days=email_config.get("days_back", 5))).strftime(
        "%d-%b-%Y"
    )
    from_address = email_config.get("from_address", DEFAULT_PUBMED_FROM)
    subject_terms = email_config.get("subject_contains_all")
    if subject_terms is None:
        subject_terms = email_config.get("subject_contains")
    if subject_terms:
        if isinstance(subject_terms, str):
            subject_terms = [subject_terms]
    else:
        subject_terms = list(DEFAULT_PUBMED_SUBJECT_TERMS)
    subject_terms = [term for term in subject_terms if term]
    subject_filters = " ".join(f'SUBJECT "{term}"' for term in subject_terms)
    criteria_parts = [f"SINCE {since_date}", f'FROM "{from_address}"']
    if subject_filters:
        criteria_parts.append(subject_filters)
    criteria = f'({" ".join(criteria_parts)})'
    status, data = client.search(None, criteria)
    if status != "OK":
        return []
    return data[0].split()


def fetch_message(client: imaplib.IMAP4_SSL, uid: bytes) -> email.message.Message:
    status, data = client.fetch(uid, "(RFC822)")
    if status != "OK":
        raise RuntimeError(f"Failed to fetch email {uid!r}")
    raw = data[0][1]
    return email.message_from_bytes(raw)


def parse_emails(client: imaplib.IMAP4_SSL, config: dict) -> tuple[list[str], str]:
    pmids = []
    search_name = ""
    for uid in search_emails(client, config):
        message = fetch_message(client, uid)
        subject = decode_mime_words(message.get("Subject", ""))
        if not search_name:
            search_name = extract_search_name(subject)
        text = extract_text_from_email(message)
        combined = f"{subject}\n{text}"
        alert_text = extract_alert_section(combined)
        pmids.extend(find_pmids(alert_text))
    return sorted(set(pmids)), search_name


def main() -> None:
    config_path = os.environ.get("CONFIG_PATH", "config.yaml")
    config = load_config(config_path)

    with connect_imap(config) as client:
        pmids, search_name = parse_emails(client, config)

    summaries = fetch_pubmed_summary(pmids)
    abstracts = fetch_pubmed_abstracts(pmids)
    client = openai_client()
    results = []
    high_if_list = config["filters"]["high_if_journals"]

    for pmid in pmids:
        article = build_article(pmid, summaries, abstracts)
        if not article.abstract:
            continue
        high_if_flag = is_high_if(article.journal, high_if_list)
        is_novel, novelty_reason = llm_novelty(client, config, article)
        if not (high_if_flag or is_novel):
            continue
        summary, strengths = llm_summary(client, config, article)
        results.append(
            ReviewResult(
                article=article,
                high_if=high_if_flag,
                is_novel=is_novel,
                novelty_reason=novelty_reason,
                summary=summary,
                strengths=strengths,
            )
        )

    sheet_name = resolve_sheet_name(config, search_name)
    append_rows(config, sheet_name, build_rows(results))


if __name__ == "__main__":
    main()
