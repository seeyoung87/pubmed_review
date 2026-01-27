import json
import logging
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Iterable
import xml.etree.ElementTree as ET

import requests
import yaml
from Bio import Entrez
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from openai import OpenAI

# Constants
DEFAULT_SHEET_NAME = "PubMed"
BATCH_SIZE = 200
SAVE_BATCH_SIZE = 10
MAX_RETRIES = 4
RETRY_DELAYS = [2, 4, 8, 16]

COLUMN_HEADERS = [
    "Date",
    "PMID",
    "Title",
    "Journal",
    "Publication Date",
    "DOI",
    "Selection Criteria",
    "Novelty Reason",
    "Summary",
    "Strengths",
]

# Default LLM Configuration
DEFAULT_LLM_CONFIG = {
    "model": "gpt-4o-mini",
    "temperature": 0.2,
    "novelty_prompt": """Determine whether the paper demonstrates genuinely novel, cutting-edge methodology.
Only mark as novel if the contribution is clearly new and meaningfully distinct from existing work.
Title: {title}
Journal: {journal}
Abstract: {abstract}""",
    "novelty_schema": {
        "name": "novelty_assessment",
        "schema": {
            "type": "object",
            "properties": {
                "is_novel": {"type": "boolean"},
                "reason": {"type": "string"}
            },
            "required": ["is_novel", "reason"],
            "additionalProperties": False
        }
    },
    "summary_prompt": """Summarize the abstract in no more than 3 lines and describe why it is strong or impactful.
Title: {title}
Journal: {journal}
Abstract: {abstract}""",
    "summary_schema": {
        "name": "summary_response",
        "schema": {
            "type": "object",
            "properties": {
                "summary": {"type": "string"},
                "strengths": {"type": "string"}
            },
            "required": ["summary", "strengths"],
            "additionalProperties": False
        }
    }
}

# Default PubMed Configuration
DEFAULT_PUBMED_CONFIG = {
    "retmax": 200,
    "datetype": "edat",
}

# Default Workflow Configuration
DEFAULT_WORKFLOW_CONFIG = {
    "fallback_days": 3,  # Used when last run time is unavailable
}

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
    """Load configuration from YAML file and merge with defaults."""
    LOGGER.info("Loading config from %s", path)
    with open(path, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    # Merge with defaults
    if "llm" not in config:
        config["llm"] = {}
    for key, value in DEFAULT_LLM_CONFIG.items():
        config["llm"].setdefault(key, value)

    if "pubmed" not in config:
        config["pubmed"] = {}
    for key, value in DEFAULT_PUBMED_CONFIG.items():
        config["pubmed"].setdefault(key, value)

    if "workflow" not in config:
        config["workflow"] = {}
    for key, value in DEFAULT_WORKFLOW_CONFIG.items():
        config["workflow"].setdefault(key, value)

    return config


def validate_config(config: dict) -> None:
    """Validate configuration has all required fields."""
    errors = []

    # Check pubmed section
    if "pubmed" not in config:
        errors.append("Missing required section: pubmed")
    else:
        pubmed = config["pubmed"]

        # Check email
        if not pubmed.get("email") and not os.environ.get("PUBMED_EMAIL"):
            errors.append("Missing required field: pubmed.email (or PUBMED_EMAIL env var)")

        # Check searches
        if "searches" not in pubmed or not pubmed["searches"]:
            errors.append("Missing required field: pubmed.searches (must have at least one search)")
        else:
            for i, search in enumerate(pubmed["searches"]):
                if not search.get("query"):
                    errors.append(f"Missing query in pubmed.searches[{i}]")

    # Check sheets section
    if "sheets" not in config:
        errors.append("Missing required section: sheets")
    else:
        sheets = config["sheets"]
        if not sheets.get("spreadsheet_id") and not os.environ.get("SPREADSHEET_ID"):
            errors.append("Missing required field: sheets.spreadsheet_id (or SPREADSHEET_ID env var)")

    # Check filters section
    if "filters" not in config:
        errors.append("Missing required section: filters")
    else:
        filters = config["filters"]
        if "high_if_journals" not in filters:
            errors.append("Missing required field: filters.high_if_journals")
        elif not isinstance(filters["high_if_journals"], list):
            errors.append("filters.high_if_journals must be a list")

    if errors:
        error_message = "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        raise ValueError(error_message)


def chunked(values: Iterable[str], size: int) -> Iterable[list[str]]:
    """Split iterable into chunks of specified size."""
    batch: list[str] = []
    for value in values:
        batch.append(value)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def calculate_reldate(last_run_time: str | None, fallback_days: int = 3) -> int:
    """Calculate reldate (days) from last successful run time.

    Args:
        last_run_time: ISO 8601 timestamp from GitHub Actions (e.g., "2024-01-15T02:00:00Z")
        fallback_days: Days to use if last_run_time is unavailable

    Returns:
        Number of days to search back, with 1 day buffer for safety
    """
    if not last_run_time:
        LOGGER.info("No last run time available, using fallback: %d days", fallback_days)
        return fallback_days

    try:
        # Parse ISO 8601 timestamp
        last_run = datetime.fromisoformat(last_run_time.replace("Z", "+00:00"))
        now = datetime.now(timezone.utc)
        delta_days = (now - last_run).days + 1  # +1 day buffer for safety

        # Minimum 1 day, maximum fallback_days * 2 (prevent excessive queries)
        reldate = max(1, min(delta_days, fallback_days * 2))
        LOGGER.info("Calculated reldate: %d days (last run: %s)", reldate, last_run_time)
        return reldate
    except (ValueError, TypeError) as exc:
        LOGGER.warning("Failed to parse last run time '%s': %s, using fallback", last_run_time, exc)
        return fallback_days


def retry_with_backoff(func, *args, max_retries=MAX_RETRIES, **kwargs):
    """Retry a function with exponential backoff."""
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except (requests.RequestException, Exception) as exc:
            if attempt == max_retries - 1:
                raise
            delay = RETRY_DELAYS[attempt] if attempt < len(RETRY_DELAYS) else RETRY_DELAYS[-1]
            LOGGER.warning("Attempt %d failed: %s. Retrying in %ds...", attempt + 1, exc, delay)
            time.sleep(delay)


def fetch_pubmed_data(pmids: list[str]) -> tuple[dict, dict[str, str]]:
    """Fetch both summaries and abstracts for PMIDs."""
    if not pmids:
        LOGGER.info("No PMIDs to fetch")
        return {}, {}

    summaries = {}
    abstracts = {}

    LOGGER.info("Fetching PubMed data for %d PMIDs", len(pmids))

    for batch in chunked(pmids, BATCH_SIZE):
        # Fetch summaries
        def fetch_summary():
            response = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "pubmed", "id": ",".join(batch), "retmode": "json"},
                timeout=30,
            )
            response.raise_for_status()
            return response.json().get("result", {})

        summaries.update(retry_with_backoff(fetch_summary))

        # Fetch abstracts
        def fetch_abstract():
            response = requests.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                params={"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
                timeout=30,
            )
            response.raise_for_status()
            return response.text

        xml_text = retry_with_backoff(fetch_abstract)
        root = ET.fromstring(xml_text)

        for article in root.findall(".//PubmedArticle"):
            pmid_element = article.find(".//MedlineCitation/PMID")
            if pmid_element is None or not pmid_element.text:
                continue

            abstract_texts = []
            for abstract in article.findall(".//AbstractText"):
                if abstract.text:
                    abstract_texts.append(abstract.text.strip())

            abstracts[pmid_element.text.strip()] = "\n".join(abstract_texts)

    return summaries, abstracts


def build_article(pmid: str, summary: dict, abstracts: dict[str, str]) -> Article:
    """Build Article object from PubMed data."""
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
    """Check if journal matches any high IF journal pattern (case-insensitive).

    Matching rules:
    - Exact match by default: "Nature" matches only "Nature"
    - Wildcard support: "Nature*" matches "Nature Medicine", "Nature Biotechnology", etc.
    - "*Radiology" matches "European Radiology", etc.
    - "The Lancet*" matches "The Lancet", "The Lancet Oncology", etc.
    """
    normalized = journal.lower()

    for pattern in high_if_list:
        pattern_lower = pattern.lower()

        # Check if pattern contains wildcard
        if "*" in pattern_lower:
            # Convert wildcard pattern to regex-like matching
            if pattern_lower.endswith("*"):
                # Prefix match: "Nature*" matches "Nature Medicine"
                prefix = pattern_lower[:-1]  # Remove trailing *
                if normalized.startswith(prefix):
                    return True
            elif pattern_lower.startswith("*"):
                # Suffix match: "*Radiology" matches "European Radiology"
                suffix = pattern_lower[1:]  # Remove leading *
                if normalized.endswith(suffix):
                    return True
            else:
                # Contains match: "The*Lancet" matches "The Lancet"
                # Split by * and check if all parts appear in order
                parts = pattern_lower.split("*")
                pos = 0
                match = True
                for part in parts:
                    idx = normalized.find(part, pos)
                    if idx == -1:
                        match = False
                        break
                    pos = idx + len(part)
                if match:
                    return True
        else:
            # Exact match (case-insensitive)
            if normalized == pattern_lower:
                return True

    return False


def openai_client() -> OpenAI:
    """Initialize OpenAI client."""
    LOGGER.info("Initializing OpenAI client")
    return OpenAI(api_key=os.environ["OPENAI_API_KEY"])


def llm_call(client: OpenAI, config: dict, prompt: str, schema: dict, max_tokens: int) -> dict:
    """Make LLM API call with structured output."""
    model = config["llm"]["model"]

    # Build API parameters
    params = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "response_format": {
            "type": "json_schema",
            "json_schema": schema,
        },
        "max_completion_tokens": max_tokens,
    }

    # Only add temperature for models that support it (not o1/o3/o4/gpt-5 reasoning models)
    if not any(x in model for x in ["o1", "o3", "o4", "gpt-5"]):
        params["temperature"] = config["llm"].get("temperature", 0.2)

    response = client.chat.completions.create(**params)

    # Log token usage
    usage = response.usage
    LOGGER.info(
        "LLM usage - prompt: %d, completion: %d, total: %d tokens",
        usage.prompt_tokens,
        usage.completion_tokens,
        usage.total_tokens,
    )

    content = response.choices[0].message.content
    try:
        return json.loads(content)
    except json.JSONDecodeError as exc:
        LOGGER.error("Failed to parse JSON. Raw content: %s", content)
        raise ValueError("Failed to parse LLM response") from exc


def llm_novelty(client: OpenAI, config: dict, article: Article) -> tuple[bool, str]:
    """Evaluate article novelty using LLM."""
    LOGGER.info("Evaluating novelty for PMID %s", article.pmid)
    prompt = config["llm"]["novelty_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    data = llm_call(client, config, prompt, config["llm"]["novelty_schema"], max_tokens=400)
    return bool(data.get("is_novel")), data.get("reason", "")


def llm_summary(client: OpenAI, config: dict, article: Article) -> tuple[str, str]:
    """Generate article summary using LLM."""
    LOGGER.info("Generating summary for PMID %s", article.pmid)
    prompt = config["llm"]["summary_prompt"].format(
        title=article.title, journal=article.journal, abstract=article.abstract
    )
    data = llm_call(client, config, prompt, config["llm"]["summary_schema"], max_tokens=500)
    return data.get("summary", ""), data.get("strengths", "")


def google_sheets_service():
    """Initialize Google Sheets service."""
    json_data = os.environ.get("GOOGLE_SERVICE_ACCOUNT_JSON")
    if not json_data:
        raise RuntimeError("Missing GOOGLE_SERVICE_ACCOUNT_JSON")

    LOGGER.info("Initializing Google Sheets service")
    info = json.loads(json_data)
    credentials = service_account.Credentials.from_service_account_info(
        info, scopes=["https://www.googleapis.com/auth/spreadsheets"]
    )
    return build("sheets", "v4", credentials=credentials)


def get_service_account_email() -> str | None:
    """Extract service account email from credentials."""
    json_data = os.environ.get("GOOGLE_SERVICE_ACCOUNT_JSON")
    if not json_data:
        return None
    try:
        info = json.loads(json_data)
        return info.get("client_email")
    except json.JSONDecodeError:
        return None


def get_existing_pmids(service, sheet_id: str, sheet_name: str) -> set[str]:
    """Get existing PMIDs from Google Sheets to avoid duplicates."""
    try:
        result = service.spreadsheets().values().get(
            spreadsheetId=sheet_id,
            range=f"{sheet_name}!B:B"  # Column B contains PMIDs
        ).execute()

        values = result.get("values", [])
        # Skip header row and extract PMIDs
        pmids = {row[0] for row in values[1:] if row and row[0] and row[0] != "PMID"}
        LOGGER.info("Found %d existing PMIDs in sheet", len(pmids))
        return pmids
    except HttpError as exc:
        if exc.resp and exc.resp.status == 404:
            LOGGER.info("Sheet %s not found, will create new", sheet_name)
            return set()
        LOGGER.warning("Could not fetch existing PMIDs: %s", exc)
        return set()


def ensure_headers(service, sheet_id: str, sheet_name: str) -> None:
    """Ensure sheet has column headers."""
    try:
        result = service.spreadsheets().values().get(
            spreadsheetId=sheet_id,
            range=f"{sheet_name}!A1:J1"
        ).execute()

        values = result.get("values", [])
        if values and values[0] == COLUMN_HEADERS:
            LOGGER.info("Headers already exist")
            return

        # Add or update headers
        service.spreadsheets().values().update(
            spreadsheetId=sheet_id,
            range=f"{sheet_name}!A1:J1",
            valueInputOption="RAW",
            body={"values": [COLUMN_HEADERS]}
        ).execute()
        LOGGER.info("Headers added to sheet")
    except HttpError:
        # If sheet doesn't exist, headers will be added with first data
        LOGGER.info("Will add headers with first data batch")


def append_rows(service, sheet_id: str, sheet_name: str, rows: list[list[str]]) -> None:
    """Append rows to Google Sheets with error handling."""
    if not rows:
        return

    LOGGER.info("Appending %d rows to sheet %s", len(rows), sheet_name)

    try:
        service.spreadsheets().values().append(
            spreadsheetId=sheet_id,
            range=f"{sheet_name}!A1",
            valueInputOption="RAW",
            insertDataOption="INSERT_ROWS",
            body={"values": rows},
        ).execute()
    except HttpError as exc:
        handle_sheets_error(exc)


def handle_sheets_error(exc: HttpError) -> None:
    """Handle Google Sheets API errors with helpful messages."""
    if not exc.resp or exc.resp.status != 403:
        raise

    # Parse error details
    reason = None
    activation_url = None
    try:
        content = exc.content.decode("utf-8") if isinstance(exc.content, (bytes, bytearray)) else exc.content
        payload = json.loads(content)
        error = payload.get("error", {})

        for detail in error.get("details", []) or []:
            if detail.get("@type") == "type.googleapis.com/google.rpc.ErrorInfo":
                reason = detail.get("reason")
                activation_url = detail.get("metadata", {}).get("activationUrl")
                break
    except (json.JSONDecodeError, TypeError, AttributeError):
        pass

    # Provide helpful error messages
    if reason == "SERVICE_DISABLED":
        url = activation_url or "https://console.developers.google.com/apis/api/sheets.googleapis.com/overview"
        raise RuntimeError(f"Google Sheets API is disabled. Enable it at {url}") from exc

    if reason == "PERMISSION_DENIED":
        email = get_service_account_email()
        hint = f"Share the spreadsheet with: {email}" if email else "Share the spreadsheet with service account"
        raise RuntimeError(f"Permission denied. {hint}") from exc

    raise


def build_row(result: ReviewResult) -> list[str]:
    """Build a single row for Google Sheets."""
    article = result.article
    selection = []
    if result.high_if:
        selection.append("High IF")
    if result.is_novel:
        selection.append("Novelty")

    return [
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


def process_article(
    article: Article,
    client: OpenAI,
    config: dict,
    high_if_list: list[str],
) -> ReviewResult | None:
    """Process a single article and return ReviewResult if it passes filters."""
    if not article.abstract:
        LOGGER.info("Skipping PMID %s - no abstract", article.pmid)
        return None

    # Check High IF first (cheaper than LLM call)
    high_if_flag = is_high_if(article.journal, high_if_list)

    # Only check novelty if not already High IF
    is_novel = False
    novelty_reason = ""
    if not high_if_flag:
        is_novel, novelty_reason = llm_novelty(client, config, article)
        if not is_novel:
            LOGGER.info("Skipping PMID %s - neither High IF nor Novel", article.pmid)
            return None
    else:
        novelty_reason = "Not evaluated (High IF journal)"

    # Generate summary for articles that passed filters
    summary, strengths = llm_summary(client, config, article)

    return ReviewResult(
        article=article,
        high_if=high_if_flag,
        is_novel=is_novel,
        novelty_reason=novelty_reason,
        summary=summary,
        strengths=strengths,
    )


def fetch_pubmed_ids(pubmed_config: dict, search_query: str, reldate: int) -> list[str]:
    """Search PubMed and return list of PMIDs."""
    reldate = int(os.environ.get("PUBMED_RELDATE", reldate))
    datetype = os.environ.get("PUBMED_DATETYPE", pubmed_config.get("datetype", "edat"))
    retmax = int(os.environ.get("PUBMED_RETMAX", pubmed_config.get("retmax", 200)))

    LOGGER.info(
        "Searching PubMed - query: %s, reldate: %s, datetype: %s, retmax: %s",
        search_query, reldate, datetype, retmax,
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
    LOGGER.info("Found %s PubMed records, returning %s IDs", count, len(id_list))

    return list(id_list)


def process_search(
    search: dict,
    config: dict,
    client: OpenAI,
    service,
    sheet_id: str,
    high_if_list: list[str],
    reldate: int,
) -> None:
    """Process a single search query."""
    search_name = search.get("name", "")
    search_query = search.get("query", "")
    sheet_name = search.get("sheet_name") or search_name or config.get("sheets", {}).get("sheet_name") or DEFAULT_SHEET_NAME

    if not search_query:
        LOGGER.warning("Skipping search with empty query for sheet %s", sheet_name)
        return

    LOGGER.info("Processing search: %s -> Sheet: %s", search_query[:50], sheet_name)

    # Fetch PMIDs
    pmids = fetch_pubmed_ids(config.get("pubmed", {}), search_query, reldate)
    if not pmids:
        LOGGER.info("No new PMIDs found")
        return

    # Get existing PMIDs to avoid duplicates
    existing_pmids = get_existing_pmids(service, sheet_id, sheet_name)
    new_pmids = [pmid for pmid in pmids if pmid not in existing_pmids]

    if not new_pmids:
        LOGGER.info("All PMIDs already processed (found %d duplicates)", len(pmids))
        return

    LOGGER.info("Processing %d new PMIDs (%d duplicates skipped)", len(new_pmids), len(pmids) - len(new_pmids))

    # Ensure headers exist
    ensure_headers(service, sheet_id, sheet_name)

    # Fetch article data
    summaries, abstracts = fetch_pubmed_data(new_pmids)

    # Process articles with batch saving
    batch = []
    total_saved = 0

    for pmid in new_pmids:
        article = build_article(pmid, summaries, abstracts)
        result = process_article(article, client, config, high_if_list)

        if result:
            batch.append(build_row(result))

            # Save in batches to prevent data loss
            if len(batch) >= SAVE_BATCH_SIZE:
                append_rows(service, sheet_id, sheet_name, batch)
                total_saved += len(batch)
                LOGGER.info("Saved batch (%d rows total so far)", total_saved)
                batch = []

    # Save remaining rows
    if batch:
        append_rows(service, sheet_id, sheet_name, batch)
        total_saved += len(batch)

    LOGGER.info("Search complete - saved %d articles to %s", total_saved, sheet_name)


def main() -> None:
    """Main entry point."""
    logging.basicConfig(
        level=os.environ.get("LOG_LEVEL", "INFO").upper(),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Load and validate configuration
    config_path = os.environ.get("CONFIG_PATH", "config.yaml")
    config = load_config(config_path)
    validate_config(config)

    # Check for dry-run mode
    dry_run = os.environ.get("DRY_RUN", "").lower() in ("true", "1", "yes")
    if dry_run:
        LOGGER.info("DRY RUN MODE - No API calls will be made")
        return

    # Setup PubMed
    pubmed_config = config.get("pubmed", {})
    email = os.environ.get("PUBMED_EMAIL") or pubmed_config.get("email", "")
    if not email:
        raise RuntimeError("Missing PUBMED_EMAIL or pubmed.email in config")
    Entrez.email = email.strip()

    # Get search queries
    searches = pubmed_config.get("searches")
    if not searches:
        search_query = pubmed_config.get("search_query")
        if not search_query:
            raise RuntimeError("Missing pubmed.search_query or pubmed.searches in config")
        searches = [{
            "name": pubmed_config.get("search_name", ""),
            "query": search_query,
            "sheet_name": pubmed_config.get("sheet_name"),
        }]

    # Setup date range - auto-calculate from last run time
    last_run_time = os.environ.get("LAST_RUN_TIME")
    fallback_days = config.get("workflow", {}).get("fallback_days", 3)
    reldate = pubmed_config.get("reldate") or calculate_reldate(last_run_time, fallback_days)

    # Initialize services
    client = openai_client()
    service = google_sheets_service()
    sheet_id = os.environ.get("SPREADSHEET_ID") or config["sheets"]["spreadsheet_id"]
    high_if_list = config["filters"]["high_if_journals"]

    # Process each search
    for search in searches:
        try:
            process_search(search, config, client, service, sheet_id, high_if_list, int(reldate))
        except Exception as exc:
            LOGGER.error("Failed to process search %s: %s", search.get("name", ""), exc)
            # Continue with next search instead of failing completely
            continue


if __name__ == "__main__":
    main()
