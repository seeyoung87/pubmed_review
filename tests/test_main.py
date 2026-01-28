"""Unit tests for pubmed_review.main module."""

import json
from unittest.mock import Mock, patch, MagicMock

import pytest

from pubmed_review.main import (
    Article,
    ReviewResult,
    chunked,
    is_selected_journal,
    build_article,
    build_row,
    load_config,
    validate_config,
)


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


class TestIsSelectedJournal:
    """Tests for is_selected_journal function."""

    def test_exact_match(self):
        """Exact match should work (case-insensitive)."""
        assert is_selected_journal("Nature", ["Nature", "Science"]) is True
        assert is_selected_journal("Science", ["Nature", "Science"]) is True

    def test_exact_match_only_without_wildcard(self):
        """Without wildcard, should only match exact journal names."""
        assert is_selected_journal("Nature Medicine", ["Nature"]) is False
        assert is_selected_journal("Skeletal Radiology", ["Radiology"]) is False
        assert is_selected_journal("The Lancet Oncology", ["Lancet"]) is False

    def test_case_insensitive(self):
        """Matching should be case-insensitive."""
        assert is_selected_journal("NATURE", ["nature"]) is True
        assert is_selected_journal("nature", ["NATURE"]) is True
        assert is_selected_journal("Nature Medicine", ["nature*"]) is True
        assert is_selected_journal("NATURE MEDICINE", ["nature*"]) is True

    def test_no_match(self):
        """Non-matching journals should return False."""
        assert is_selected_journal("Random Journal", ["Nature", "Science"]) is False

    def test_empty_list(self):
        """Empty list should return False."""
        assert is_selected_journal("Nature", []) is False

    def test_wildcard_prefix_match(self):
        """Wildcard at end should match journals starting with pattern."""
        assert is_selected_journal("Nature Medicine", ["Nature*"]) is True
        assert is_selected_journal("Nature Biotechnology", ["Nature*"]) is True
        assert is_selected_journal("Nature", ["Nature*"]) is True
        assert is_selected_journal("The Lancet", ["The Lancet*"]) is True
        assert is_selected_journal("The Lancet Oncology", ["The Lancet*"]) is True

    def test_wildcard_suffix_match(self):
        """Wildcard at start should match journals ending with pattern."""
        assert is_selected_journal("European Radiology", ["*Radiology"]) is True
        assert is_selected_journal("Skeletal Radiology", ["*Radiology"]) is True
        assert is_selected_journal("Radiology", ["*Radiology"]) is True

    def test_wildcard_no_false_positives(self):
        """Wildcard should not match unrelated journals."""
        assert is_selected_journal("Science", ["Nature*"]) is False
        assert is_selected_journal("Radiology AI", ["*Radiology"]) is False
        assert is_selected_journal("Neuroradiology", ["Radiology*"]) is False

    def test_multiple_patterns(self):
        """Should match if any pattern matches."""
        patterns = ["Nature", "Science*", "*Radiology"]
        assert is_selected_journal("Nature", patterns) is True
        assert is_selected_journal("Science Advances", patterns) is True
        assert is_selected_journal("European Radiology", patterns) is True
        assert is_selected_journal("Random Journal", patterns) is False


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

    def test_selected_journal_and_novel(self):
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
            selected_journal=True,
            is_novel=True,
            summary="Summary text",
            comment="Noteworthy finding"
        )

        row = build_row(result)

        assert len(row) == 6
        assert row[1] == "Test"    # Title
        assert row[2] == "Nature"  # Journal
        assert row[3] == "Selected Journal, Novelty"  # Selection Criteria
        assert row[4] == "Summary text"
        assert row[5] == "Noteworthy finding"

    def test_selected_journal_only(self):
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
            selected_journal=True,
            is_novel=False,
            summary="Summary",
            comment="Comment"
        )

        row = build_row(result)
        assert row[3] == "Selected Journal"

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
            selected_journal=False,
            is_novel=True,
            summary="Summary",
            comment="Comment"
        )

        row = build_row(result)
        assert row[3] == "Novelty"


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
        assert config["llm"]["model"] == "gpt-5-mini"

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
  model: "gpt-5"
""")

        config = load_config(str(config_file))

        # Check overrides work
        assert config["llm"]["model"] == "gpt-5"
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
