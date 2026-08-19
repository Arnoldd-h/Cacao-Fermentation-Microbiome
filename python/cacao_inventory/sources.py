"""Read-only clients for authoritative NCBI and ENA metadata services."""

from __future__ import annotations

import csv
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any


ENA_FIELDS = (
    "study_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "run_accession",
    "experiment_title",
    "sample_alias",
    "sample_title",
    "sample_description",
    "scientific_name",
    "tax_id",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "base_count",
    "read_count",
    "fastq_bytes",
    "fastq_ftp",
    "fastq_md5",
    "country",
    "collection_date",
    "cultivar",
    "isolation_source",
    "first_public",
    "last_updated",
)


class MetadataSourceError(RuntimeError):
    """Raised when an authoritative metadata source cannot be used safely."""


class MetadataClient:
    """Fetch public project summaries and run reports without FASTQ downloads."""

    def __init__(self, inventory_config: dict[str, Any]) -> None:
        self.ncbi_search_url = str(inventory_config["ncbi_esearch_url"])
        self.ncbi_url = str(inventory_config["ncbi_esummary_url"])
        self.ena_url = str(inventory_config["ena_filereport_url"])
        self.timeout = int(inventory_config.get("http_timeout_seconds", 60))
        self.retries = int(inventory_config.get("http_retries", 3))
        self.user_agent = str(inventory_config.get("user_agent", "cacao-inventory/0.1"))
        if self.retries < 1:
            raise ValueError("http_retries must be at least one")

    def _get_text(self, base_url: str, params: dict[str, str]) -> str:
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        request = urllib.request.Request(url, headers={"User-Agent": self.user_agent})
        last_error: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read().decode("utf-8")
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(0.5 * attempt)
        raise MetadataSourceError(f"Metadata request failed after {self.retries} attempts: {url}") from last_error

    def _search_uids(self, term: str, retmax: int) -> list[str]:
        search_text = self._get_text(
            self.ncbi_search_url,
            {
                "db": "bioproject",
                "term": term,
                "retmode": "json",
                "retmax": str(retmax),
            },
        )
        try:
            search_payload = json.loads(search_text)
            uids = search_payload["esearchresult"]["idlist"]
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            raise MetadataSourceError("NCBI BioProject response is not valid ESearch JSON") from exc
        if not uids:
            raise MetadataSourceError("NCBI ESearch returned no BioProject UIDs")
        return [str(uid) for uid in uids]

    def _fetch_summaries_by_uid(self, uids: list[str]) -> dict[str, dict[str, Any]]:
        text = self._get_text(
            self.ncbi_url,
            {"db": "bioproject", "id": ",".join(str(uid) for uid in uids), "retmode": "json"},
        )
        try:
            payload = json.loads(text)
            result = payload["result"]
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            raise MetadataSourceError("NCBI BioProject response is not valid ESummary JSON") from exc

        summaries: dict[str, dict[str, Any]] = {}
        for uid in result.get("uids", []):
            record = result.get(str(uid), {})
            accession = str(record.get("project_acc", ""))
            if accession:
                summaries[accession] = record
        return summaries

    def fetch_bioproject_summaries(self, accessions: list[str]) -> dict[str, dict[str, Any]]:
        """Return NCBI BioProject summaries keyed by accession."""

        uids = self._search_uids(" OR ".join(accessions), len(accessions))
        summaries = self._fetch_summaries_by_uid(uids)
        missing = sorted(set(accessions) - set(summaries))
        if missing:
            raise MetadataSourceError(f"NCBI returned no valid BioProject summary for: {', '.join(missing)}")
        return summaries

    def search_bioproject_summaries(self, term: str, retmax: int = 200) -> dict[str, dict[str, Any]]:
        """Search NCBI BioProject and return all summaries found for a documented term."""

        uids = self._search_uids(term, retmax)
        summaries = self._fetch_summaries_by_uid(uids)
        if len(summaries) != len(uids):
            raise MetadataSourceError(
                f"NCBI returned {len(summaries)} summaries for {len(uids)} discovered UIDs"
            )
        return summaries

    def fetch_ena_runs(self, accession: str) -> list[dict[str, str]]:
        """Return all ENA read-run rows for a BioProject accession."""

        text = self._get_text(
            self.ena_url,
            {
                "accession": accession,
                "result": "read_run",
                "fields": ",".join(ENA_FIELDS),
                "format": "tsv",
                "download": "true",
            },
        )
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        headers = set(reader.fieldnames or [])
        missing = sorted(set(ENA_FIELDS) - headers)
        if missing:
            raise MetadataSourceError(
                f"ENA schema changed or is incomplete for {accession}; missing: {', '.join(missing)}"
            )
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        for row in rows:
            if not row.get("run_accession"):
                raise MetadataSourceError(f"ENA returned a row without run_accession for {accession}")
        return rows
