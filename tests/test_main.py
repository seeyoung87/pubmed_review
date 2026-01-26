"""Unit tests for pubmed_review.main module."""

import json
from unittest.mock import Mock, patch, MagicMock

import pytest

from pubmed_review.main import (
    Article,
    ReviewResult,
    chunked,
    is_high_if,
    build_article,
    build_row,
    load_config,
    validate_config,
    calculate_reldate,
)


class TestCalculateReldate:
    """Tests for calculate_reldate function."""

    def test_no_last_run_time(self):
        """Should return fallback_days when no last run time."""
        assert calculate_reldate(None, fallback_days=3) == 3
        assert calculate_reldate("", fallback_days=5) == 5

    def test_recent_run(self):
        """Should calculate days from last run with buffer."""
        from datetime import datetime, timezone, timedelta

        # 1 day ago
        one_day_ago = (datetime.now(timezone.utc) - timedelta(days=1)).isoformat()
        result = calculate_reldate(one_day_ago, fallback_days=3)
        assert result >= 1  # At least 1 day
        assert result <= 3  # Within expected range

    def test_max_cap(self):
        """Should cap at fallback_days * 2."""
        from datetime import datetime, timezone, timedelta

        # 30 days ago
        thirty_days_ago = (datetime.now(timezone.utc) - timedelta(days=30)).isoformat()
        result = calculate_reldate(thirty_days_ago, fallback_days=3)
        assert result == 6  # Capped at 3 * 2

    def test_invalid_timestamp(self):
        """Should return fallback_days for invalid timestamp."""
        assert calculate_reldate("invalid", fallback_days=3) == 3
        assert calculate_reldate("2024-13-45", fallback_days=5) == 5

    def test_github_format(self):
        """Should handle GitHub API timestamp format."""
        from datetime import datetime, timezone, timedelta

        # GitHub format with Z suffix
        two_days_ago = (datetime.now(timezone.utc) - timedelta(days=2)).strftime("%Y-%m-%dT%H:%M:%SZ")
        result = calculate_reldate(two_days_ago, fallback_days=7)
        assert result >= 2  # At least 2 days + buffer
        assert result <= 7  # Within fallback


class TestChunked:
    """Tests for chunked function."""

    def test_empty_list(self):
        result = list(chunked([], 3))
        assert result == []

    def test_exact_chunks(self):
        result = list(chunked(["a", "b", "c", "d", "e", "f"], 3))
        assert result == [["a", "b", "c"], ["d", "e", "f"]]

    def test_partial_chunk(self):
        result = list(chunked(["a", "b", "c", "d", "e"], 3))
        assert result == [["a", "b", "c"], ["d", "e"]]

    def test_single_chunk(self):
        result = list(chunked(["a", "b"], 5))
        assert result == [["a", "b"]]


class TestIsHighIF:
    """Tests for is_high_if function."""

    def test_exact_match(self):
        """Exact match should work (case-insensitive)."""
        assert is_high_if("Nature", ["Nature", "Science"]) is True
        assert is_high_if("Science", ["Nature", "Science"]) is True

    def test_exact_match_only_without_wildcard(self):
        """Without wildcard, should only match exact journal names."""
        assert is_high_if("Nature Medicine", ["Nature"]) is False
        assert is_high_if("Skeletal Radiology", ["Radiology"]) is False
        assert is_high_if("The Lancet Oncology", ["Lancet"]) is False

    def test_case_insensitive(self):
        """Matching should be case-insensitive."""
        assert is_high_if("NATURE", ["nature"]) is True
        assert is_high_if("nature", ["NATURE"]) is True
        assert is_high_if("Nature Medicine", ["nature*"]) is True
        assert is_high_if("NATURE MEDICINE", ["nature*"]) is True

    def test_no_match(self):
        """Non-matching journals should return False."""
        assert is_high_if("Random Journal", ["Nature", "Science"]) is False

    def test_empty_list(self):
        """Empty list should return False."""
        assert is_high_if("Nature", []) is False

    def test_wildcard_prefix_match(self):
        """Wildcard at end should match journals starting with pattern."""
        assert is_high_if("Nature Medicine", ["Nature*"]) is True
        assert is_high_if("Nature Biotechnology", ["Nature*"]) is True
        assert is_high_if("Nature", ["Nature*"]) is True
        assert is_high_if("The Lancet", ["The Lancet*"]) is True
        assert is_high_if("The Lancet Oncology", ["The Lancet*"]) is True

    def test_wildcard_suffix_match(self):
        """Wildcard at start should match journals ending with pattern."""
        assert is_high_if("European Radiology", ["*Radiology"]) is True
        assert is_high_if("Skeletal Radiology", ["*Radiology"]) is True
        assert is_high_if("Radiology", ["*Radiology"]) is True

    def test_wildcard_no_false_positives(self):
        """Wildcard should not match unrelated journals."""
        assert is_high_if("Science", ["Nature*"]) is False
        assert is_high_if("Radiology AI", ["*Radiology"]) is False
        assert is_high_if("Neuroradiology", ["Radiology*"]) is False

    def test_multiple_patterns(self):
        """Should match if any pattern matches."""
        patterns = ["Nature", "Science*", "*Radiology"]
        assert is_high_if("Nature", patterns) is True
        assert is_high_if("Science Advances", patterns) is True
        assert is_high_if("European Radiology", patterns) is True
        assert is_high_if("Random Journal", patterns) is False


class TestBuildArticle:
    """Tests for build_article function."""

    def test_complete_article(self):
        pmid = "12345"
        summary = {
            "12345": {
                "title": "Test Article",
                "fulljournalname": "Nature",
                "pubdate": "2026 Jan",
                "articleids": [
                    {"idtype": "doi", "value": "10.1038/test"},
                    {"idtype": "pmid", "value": "12345"}
                ],
                "authors": [
                    {"name": "Smith J"},
                    {"name": "Doe A"}
                ]
            }
        }
        abstracts = {"12345": "This is a test abstract."}

        article = build_article(pmid, summary, abstracts)

        assert article.pmid == "12345"
        assert article.title == "Test Article"
        assert article.journal == "Nature"
        assert article.pub_date == "2026 Jan"
        assert article.doi == "10.1038/test"
        assert article.authors == "Smith J; Doe A"
        assert article.abstract == "This is a test abstract."
        assert article.url == "https://pubmed.ncbi.nlm.nih.gov/12345/"

    def test_missing_fields(self):
        pmid = "12345"
        summary = {"12345": {}}
        abstracts = {}

        article = build_article(pmid, summary, abstracts)

        assert article.pmid == "12345"
        assert article.title == ""
        assert article.journal == ""
        assert article.doi == ""
        assert article.abstract == ""

    def test_missing_pmid_in_summary(self):
        article = build_article("99999", {}, {})
        assert article.pmid == "99999"
        assert article.title == ""


class TestBuildRow:
    """Tests for build_row function."""

    def test_high_if_and_novel(self):
        article = Article(
            pmid="12345",
            title="Test",
            journal="Nature",
            pub_date="2026 Jan",
            doi="10.1038/test",
            authors="Smith J",
            abstract="Abstract",
            url="https://pubmed.ncbi.nlm.nih.gov/12345/"
        )
        result = ReviewResult(
            article=article,
            high_if=True,
            is_novel=True,
            novelty_reason="Novel method",
            summary="Summary text",
            strengths="Strong validation"
        )

        row = build_row(result)

        assert row[1] == "12345"  # PMID
        assert row[2] == "Test"   # Title
        assert row[3] == "Nature" # Journal
        assert row[6] == "High IF, Novelty"  # Selection Criteria
        assert row[7] == "Novel method"
        assert row[8] == "Summary text"
        assert row[9] == "Strong validation"

    def test_high_if_only(self):
        article = Article(
            pmid="12345",
            title="Test",
            journal="Nature",
            pub_date="2026 Jan",
            doi="10.1038/test",
            authors="Smith J",
            abstract="Abstract",
            url="https://pubmed.ncbi.nlm.nih.gov/12345/"
        )
        result = ReviewResult(
            article=article,
            high_if=True,
            is_novel=False,
            novelty_reason="Not evaluated",
            summary="Summary",
            strengths="Strengths"
        )

        row = build_row(result)
        assert row[6] == "High IF"

    def test_novel_only(self):
        article = Article(
            pmid="12345",
            title="Test",
            journal="Other Journal",
            pub_date="2026 Jan",
            doi="10.1038/test",
            authors="Smith J",
            abstract="Abstract",
            url="https://pubmed.ncbi.nlm.nih.gov/12345/"
        )
        result = ReviewResult(
            article=article,
            high_if=False,
            is_novel=True,
            novelty_reason="Novel approach",
            summary="Summary",
            strengths="Strengths"
        )

        row = build_row(result)
        assert row[6] == "Novelty"


class TestLoadConfig:
    """Tests for load_config function."""

    def test_load_minimal_config(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
pubmed:
  email: "test@example.com"
  searches:
    - query: "test query"
      sheet_name: "Test"
sheets:
  spreadsheet_id: "test123"
filters:
  high_if_journals: ["Nature"]
""")

        config = load_config(str(config_file))

        # Check defaults are merged
        assert "llm" in config
        assert config["llm"]["model"] == "gpt-4o-mini"
        assert config["llm"]["temperature"] == 0.2
        assert "workflow" in config
        assert config["workflow"]["fallback_days"] == 3

    def test_load_with_overrides(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("""
pubmed:
  email: "test@example.com"
  retmax: 50
  searches:
    - query: "test"
      sheet_name: "Test"
sheets:
  spreadsheet_id: "test123"
filters:
  high_if_journals: ["Nature"]
llm:
  model: "gpt-4"
  temperature: 0.5
workflow:
  fallback_days: 7
""")

        config = load_config(str(config_file))

        # Check overrides work
        assert config["llm"]["model"] == "gpt-4"
        assert config["llm"]["temperature"] == 0.5
        assert config["workflow"]["fallback_days"] == 7
        assert config["pubmed"]["retmax"] == 50


class TestValidateConfig:
    """Tests for validate_config function."""

    def test_valid_config(self):
        config = {
            "pubmed": {
                "email": "test@example.com",
                "searches": [
                    {"query": "test", "sheet_name": "Test"}
                ]
            },
            "sheets": {
                "spreadsheet_id": "test123"
            },
            "filters": {
                "high_if_journals": ["Nature"]
            }
        }
        # Should not raise
        validate_config(config)

    def test_missing_email(self):
        config = {
            "pubmed": {
                "searches": [{"query": "test", "sheet_name": "Test"}]
            },
            "sheets": {"spreadsheet_id": "test123"},
            "filters": {"high_if_journals": ["Nature"]}
        }
        with pytest.raises(ValueError, match="pubmed.email"):
            validate_config(config)

    def test_missing_searches(self):
        config = {
            "pubmed": {
                "email": "test@example.com"
            },
            "sheets": {"spreadsheet_id": "test123"},
            "filters": {"high_if_journals": ["Nature"]}
        }
        with pytest.raises(ValueError, match="pubmed.searches"):
            validate_config(config)

    def test_empty_searches(self):
        config = {
            "pubmed": {
                "email": "test@example.com",
                "searches": []
            },
            "sheets": {"spreadsheet_id": "test123"},
            "filters": {"high_if_journals": ["Nature"]}
        }
        with pytest.raises(ValueError, match="at least one search"):
            validate_config(config)

    def test_missing_spreadsheet_id(self):
        config = {
            "pubmed": {
                "email": "test@example.com",
                "searches": [{"query": "test", "sheet_name": "Test"}]
            },
            "sheets": {},
            "filters": {"high_if_journals": ["Nature"]}
        }
        # Should not raise if SPREADSHEET_ID env var could be set
        # But in test, it should raise
        with patch.dict('os.environ', {}, clear=True):
            with pytest.raises(ValueError, match="spreadsheet_id"):
                validate_config(config)

    def test_missing_high_if_journals(self):
        config = {
            "pubmed": {
                "email": "test@example.com",
                "searches": [{"query": "test", "sheet_name": "Test"}]
            },
            "sheets": {"spreadsheet_id": "test123"},
            "filters": {}
        }
        with pytest.raises(ValueError, match="high_if_journals"):
            validate_config(config)

    def test_search_missing_query(self):
        config = {
            "pubmed": {
                "email": "test@example.com",
                "searches": [
                    {"sheet_name": "Test"}
                ]
            },
            "sheets": {"spreadsheet_id": "test123"},
            "filters": {"high_if_journals": ["Nature"]}
        }
        with pytest.raises(ValueError, match="query"):
            validate_config(config)
