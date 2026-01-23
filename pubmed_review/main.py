import json
import logging
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Iterable
import xml.etree.ElementTree as ET

import requests
import yaml
from Bio import Entrez
from google.oauth2 import service_account
from googleapiclient.discovery import build
from openai import OpenAI

DEFAULT_SHEET_NAME = "PubMed"

LOGGER = logging.getLogger(__name__)


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
    LOGGER.info("Loading config from %s", path)
    with open(path, "r", encoding="utf-8") as file:
        return yaml.safe_load(file)


def chunked(values: Iterable[str], size: int) -> Iterable[list[str]]:
    batch: list[str] = []
    for value in values:
        batch.append(value)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def fetch_pubmed_summary(pmids: Iterable[str], batch_size: int = 200) -> dict:
    pmid_list = list(pmids)
    if not pmid_list:
        LOGGER.info("No PMIDs found. Skipping PubMed summary fetch.")
        return {}
    summaries: dict = {}
    LOGGER.info("Fetching PubMed summary for %d PMIDs", len(pmid_list))
    for batch in chunked(pmid_list, batch_size):
        response = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "pubmed", "id": ",".join(batch), "retmode": "json"},
            timeout=30,
        )
        response.raise_for_status()
        summaries.update(response.json().get("result", {}))
    return summaries


def fetch_pubmed_abstracts(pmids: Iterable[str], batch_size: int = 200) -> dict[str, str]:
    pmid_list = list(pmids)
    if not pmid_list:
        LOGGER.info("No PMIDs found. Skipping PubMed abstract fetch.")
        return {}
    abstracts_by_pmid: dict[str, str] = {}
    LOGGER.info("Fetching PubMed abstracts for %d PMIDs", len(pmid_list))
    for batch in chunked(pmid_list, batch_size):
        response = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
            timeout=30,
        )
        response.raise_for_status()
        root = ET.fromstring(response.text)
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
    LOGGER.info("Initializing OpenAI client")
    return OpenAI(api_key=os.environ["OPENAI_API_KEY"])


def normalize_env_value(value: str) -> str:
    return value.replace("\u00a0", " ").strip()


def parse_json_response(content: str, context: str) -> dict:
    try:
        return json.loads(content)
    except json.JSONDecodeError as exc:
        LOGGER.error("Failed to parse JSON for %s. Raw content: %s", context, content)
        raise ValueError(f"Failed to parse JSON for {context}") from exc


def llm_novelty(client: OpenAI, config: dict, article: Article) -> tuple[bool, str]:
    LOGGER.info("Evaluating novelty for PMID %s", article.pmid)
    prompt = config["llm"]["novelty_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    novelty_schema = config["llm"]["novelty_schema"]
    response = client.chat.completions.create(
        model=config["llm"]["model"],
        temperature=config["llm"].get("temperature", 0.2),
        messages=[{"role": "user", "content": prompt}],
        response_format={
            "type": "json_schema",
            "json_schema": novelty_schema,
        },
        max_tokens=400,
    )
    content = response.choices[0].message.content
    data = parse_json_response(content, "novelty response")
    return bool(data.get("is_novel")), data.get("reason", "")


def llm_summary(client: OpenAI, config: dict, article: Article) -> tuple[str, str]:
    LOGGER.info("Generating summary for PMID %s", article.pmid)
    prompt = config["llm"]["summary_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    summary_schema = config["llm"]["summary_schema"]
    response = client.chat.completions.create(
        model=config["llm"]["model"],
        temperature=config["llm"].get("temperature", 0.2),
        messages=[{"role": "user", "content": prompt}],
        response_format={
            "type": "json_schema",
            "json_schema": summary_schema,
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
        LOGGER.info("Using GOOGLE_SERVICE_ACCOUNT_JSON for Sheets auth")
        info = json.loads(json_data)
    elif file_path:
        LOGGER.info("Using GOOGLE_SERVICE_ACCOUNT_FILE for Sheets auth")
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
        LOGGER.info("Using search name for sheet: %s", search_name)
        return search_name
    sheet_config = config.get("sheets", {})
    resolved = sheet_config.get("sheet_name") or DEFAULT_SHEET_NAME
    LOGGER.info("Using default sheet name: %s", resolved)
    return resolved


def append_rows(config: dict, sheet_name: str, rows: list[list[str]]) -> None:
    if not rows:
        LOGGER.info("No rows to append to Sheets. Skipping update.")
        return
    service = google_sheets_service()
    sheet_id = os.environ.get("SPREADSHEET_ID") or config["sheets"]["spreadsheet_id"]
    LOGGER.info("Appending %d rows to sheet %s", len(rows), sheet_name)
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


def load_pubmed_settings(config: dict) -> tuple[dict, str]:
    pubmed_config = config.get("pubmed", {})
    email_address = os.environ.get("PUBMED_EMAIL") or pubmed_config.get("email")
    if not email_address:
        raise RuntimeError("Missing PUBMED_EMAIL or pubmed.email in config")
    search_query = os.environ.get("PUBMED_SEARCH_QUERY") or pubmed_config.get("search_query")
    if not search_query:
        raise RuntimeError("Missing PUBMED_SEARCH_QUERY or pubmed.search_query in config")
    return pubmed_config, search_query


def fetch_pubmed_ids(pubmed_config: dict, search_query: str, reldate: int) -> list[str]:
    email_address = os.environ.get("PUBMED_EMAIL") or pubmed_config.get("email", "")
    Entrez.email = normalize_env_value(email_address)
    reldate = int(os.environ.get("PUBMED_RELDATE", reldate))
    datetype = os.environ.get("PUBMED_DATETYPE", pubmed_config.get("datetype", "edat"))
    retmax = int(os.environ.get("PUBMED_RETMAX", pubmed_config.get("retmax", 200)))
    LOGGER.info(
        "Searching PubMed with query=%s reldate=%s datetype=%s retmax=%s",
        search_query,
        reldate,
        datetype,
        retmax,
    )
    with Entrez.esearch(
        db="pubmed",
        term=search_query,
        reldate=reldate,
        datetype=datetype,
        retmax=retmax,
    ) as handle:
        record = Entrez.read(handle)
    count = record.get("Count", "0")
    id_list = record.get("IdList", [])
    LOGGER.info("Found %s PubMed records. IDs returned: %s", count, len(id_list))
    return list(id_list)


def main() -> None:
    logging.basicConfig(
        level=os.environ.get("LOG_LEVEL", "INFO").upper(),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    config_path = os.environ.get("CONFIG_PATH", "config.yaml")
    config = load_config(config_path)

    pubmed_config, search_query = load_pubmed_settings(config)
    schedule_days = config.get("workflow", {}).get("schedule_days", 1)
    reldate = pubmed_config.get("reldate")
    if reldate is None:
        reldate = schedule_days
    pmids = fetch_pubmed_ids(pubmed_config, search_query, int(reldate))
    search_name = pubmed_config.get("search_name", "")

    if not pmids:
        LOGGER.info("No PubMed records found. Exiting without updates.")
        return

    summaries = fetch_pubmed_summary(pmids)
    abstracts = fetch_pubmed_abstracts(pmids)
    client = openai_client()
    results = []
    high_if_list = config["filters"]["high_if_journals"]

    for pmid in pmids:
        article = build_article(pmid, summaries, abstracts)
        if not article.abstract:
            LOGGER.info("Skipping PMID %s because abstract is missing", pmid)
            continue
        high_if_flag = is_high_if(article.journal, high_if_list)
        is_novel, novelty_reason = llm_novelty(client, config, article)
        if not (high_if_flag or is_novel):
            LOGGER.info("Skipping PMID %s because it is neither High IF nor Novel", pmid)
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
