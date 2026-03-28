"""
pubmed_fetch.py  —  Phase 0, Step 1
Queries PubMed (NCBI E-utilities) and EuropePMC (REST API) for research papers
on Pseudomonas aeruginosa redox biology and antimicrobial resistance.

Target: >= 1000 unique papers (deduplicated by PMID/DOI).

Output:
    literature_mining/raw_papers/pubmed_raw.json
    literature_mining/raw_papers/epmc_raw.json

Usage:
    python pubmed_fetch.py
    python pubmed_fetch.py --max-per-query 200
"""

import argparse
import json
import time
import requests
from pathlib import Path

# ---------------------------------------------------------------------------
# Search queries — designed to cover all redox/AMR angles in P. aeruginosa
# Each query targets a distinct biological mechanism to maximise coverage
# and minimise redundancy between queries.
# ---------------------------------------------------------------------------
SEARCH_QUERIES = [
    # Core: redox biology in P. aeruginosa
    "Pseudomonas aeruginosa redox",
    "Pseudomonas aeruginosa oxidative stress",
    "Pseudomonas aeruginosa reactive oxygen species",
    "Pseudomonas aeruginosa ROS resistance",

    # Antimicrobial resistance — broad
    "Pseudomonas aeruginosa antimicrobial resistance mechanism",
    "Pseudomonas aeruginosa antibiotic resistance oxidative",
    "Pseudomonas aeruginosa multidrug resistance redox",

    # Redox enzymes and detoxification systems
    "Pseudomonas aeruginosa catalase peroxidase",
    "Pseudomonas aeruginosa superoxide dismutase",
    "Pseudomonas aeruginosa glutathione",
    "Pseudomonas aeruginosa thioredoxin",
    "Pseudomonas aeruginosa glutaredoxin",
    "Pseudomonas aeruginosa ferredoxin",
    "Pseudomonas aeruginosa alkyl hydroperoxide reductase",
    "Pseudomonas aeruginosa OxyR",

    # Electron transport and energy metabolism
    "Pseudomonas aeruginosa electron transport chain",
    "Pseudomonas aeruginosa NADH dehydrogenase",
    "Pseudomonas aeruginosa cytochrome oxidase",
    "Pseudomonas aeruginosa proton motive force antibiotic",
    "Pseudomonas aeruginosa terminal oxidase",

    # Efflux pumps — redox/PMF-driven
    "Pseudomonas aeruginosa MexAB-OprM efflux",
    "Pseudomonas aeruginosa efflux pump resistance",
    "Pseudomonas aeruginosa RND efflux oxidative",

    # Biofilm — redox microenvironment
    "Pseudomonas aeruginosa biofilm redox",
    "Pseudomonas aeruginosa biofilm electron acceptor",
    "Pseudomonas aeruginosa phenazine electron shuttle",
    "Pseudomonas aeruginosa pyocyanin redox",

    # Specific resistance-linked redox proteins
    "Pseudomonas aeruginosa KatA catalase",
    "Pseudomonas aeruginosa KatB catalase",
    "Pseudomonas aeruginosa AhpC",
    "Pseudomonas aeruginosa MvfR redox",
    "Pseudomonas aeruginosa SoxR",
    "Pseudomonas aeruginosa Fur iron redox",

    # Iron/metal redox — critical in P. aeruginosa virulence
    "Pseudomonas aeruginosa iron acquisition redox",
    "Pseudomonas aeruginosa siderophore oxidative",
    "Pseudomonas aeruginosa Fenton reaction",

    # Quorum sensing + redox
    "Pseudomonas aeruginosa quorum sensing oxidative stress",
    "Pseudomonas aeruginosa las rhl redox",

    # Broader AMR + redox synergy
    "Pseudomonas aeruginosa beta-lactam resistance oxidative",
    "Pseudomonas aeruginosa aminoglycoside resistance redox",
    "Pseudomonas aeruginosa fluoroquinolone oxidative stress",
    "Pseudomonas aeruginosa persister cell redox",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
NCBI_RATE_LIMIT = 0.34   # seconds between requests (3 req/s max allowed by NCBI)
OUTPUT_DIR = Path(__file__).parent / "raw_papers"


# ---------------------------------------------------------------------------
# PubMed fetch
# ---------------------------------------------------------------------------
def fetch_pubmed(query: str, max_results: int) -> list:
    """Fetch papers from PubMed for a single query. Returns list of records."""
    records = []

    # Step 1: esearch — get PMIDs
    search_params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "usehistory": "y",
    }
    try:
        resp = requests.get(
            f"{NCBI_BASE}/esearch.fcgi", params=search_params, timeout=30
        )
        resp.raise_for_status()
        search_data = resp.json()
    except Exception as e:
        print(f"  [PubMed] esearch failed for '{query}': {e}")
        return records

    time.sleep(NCBI_RATE_LIMIT)

    pmids = search_data.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return records

    web_env = search_data["esearchresult"].get("webenv", "")
    query_key = search_data["esearchresult"].get("querykey", "")

    # Step 2: esummary — get metadata for all PMIDs in one call
    fetch_params = {
        "db": "pubmed",
        "query_key": query_key,
        "WebEnv": web_env,
        "retmax": max_results,
        "retmode": "json",
    }
    try:
        resp = requests.get(
            f"{NCBI_BASE}/esummary.fcgi", params=fetch_params, timeout=30
        )
        resp.raise_for_status()
        fetch_data = resp.json()
    except Exception as e:
        print(f"  [PubMed] esummary failed for '{query}': {e}")
        return records

    time.sleep(NCBI_RATE_LIMIT)

    result = fetch_data.get("result", {})
    for pmid in pmids:
        paper = result.get(pmid, {})
        if not paper:
            continue
        records.append({
            "pmid": pmid,
            "doi": next(
                (
                    id_obj.get("value", "")
                    for id_obj in paper.get("articleids", [])
                    if id_obj.get("idtype") == "doi"
                ),
                "",
            ),
            "title": paper.get("title", ""),
            "authors": [a.get("name", "") for a in paper.get("authors", [])],
            "journal": paper.get("fulljournalname", ""),
            "year": paper.get("pubdate", "")[:4],
            "abstract": "",   # esummary doesn't include abstracts
            "source": "pubmed",
            "query": query,
        })

    return records


# ---------------------------------------------------------------------------
# EuropePMC fetch
# ---------------------------------------------------------------------------
def fetch_epmc(query: str, max_results: int) -> list:
    """Fetch papers from EuropePMC for a single query. Returns list of records."""
    records = []
    page_size = min(max_results, 100)
    fetched = 0
    cursor = "*"

    while fetched < max_results:
        params = {
            "query": query,
            "format": "json",
            "pageSize": page_size,
            "cursorMark": cursor,
            "resultType": "core",
        }
        try:
            resp = requests.get(EPMC_BASE, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"  [EuropePMC] fetch failed for '{query}': {e}")
            break

        results = data.get("resultList", {}).get("result", [])
        if not results:
            break

        for paper in results:
            records.append({
                "pmid": paper.get("pmid", ""),
                "doi": paper.get("doi", ""),
                "title": paper.get("title", ""),
                "authors": [
                    a.get("fullName", "")
                    for a in paper.get("authorList", {}).get("author", [])
                ],
                "journal": paper.get("journalTitle", ""),
                "year": str(paper.get("pubYear", "")),
                "abstract": paper.get("abstractText", ""),
                "source": "europepmc",
                "query": query,
            })

        fetched += len(results)
        next_cursor = data.get("nextCursorMark", "")
        if not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        time.sleep(0.2)

    return records


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------
def deduplicate(records: list) -> list:
    """
    Deduplicate by PMID first, then by DOI.
    When a duplicate is found, keep the record with an abstract (EuropePMC)
    over one without (PubMed esummary).
    """
    seen_pmids = {}
    seen_dois = {}
    unique = []

    def prefer_new(new, existing):
        return bool(new.get("abstract")) and not bool(existing.get("abstract"))

    for rec in records:
        pmid = rec.get("pmid", "").strip()
        doi = rec.get("doi", "").strip().lower()

        if pmid and pmid in seen_pmids:
            if prefer_new(rec, seen_pmids[pmid]):
                idx = unique.index(seen_pmids[pmid])
                unique[idx] = rec
                seen_pmids[pmid] = rec
            continue

        if doi and doi in seen_dois:
            if prefer_new(rec, seen_dois[doi]):
                idx = unique.index(seen_dois[doi])
                unique[idx] = rec
                seen_dois[doi] = rec
            continue

        unique.append(rec)
        if pmid:
            seen_pmids[pmid] = rec
        if doi:
            seen_dois[doi] = rec

    return unique


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch P. aeruginosa redox/AMR papers from PubMed and EuropePMC. "
            "Target: >= 1000 unique papers."
        )
    )
    parser.add_argument(
        "--max-per-query",
        type=int,
        default=150,
        help="Maximum results per query per source (default: 150).",
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pubmed_all = []
    epmc_all = []
    total_queries = len(SEARCH_QUERIES)

    print(f"Starting fetch: {total_queries} queries x 2 sources, "
          f"max {args.max_per_query} per query.")
    print(f"Theoretical max before dedup: "
          f"{total_queries * args.max_per_query * 2:,} records.\n")

    for i, query in enumerate(SEARCH_QUERIES, 1):
        print(f"[{i:02d}/{total_queries}] \"{query}\"")

        pm = fetch_pubmed(query, args.max_per_query)
        print(f"         PubMed:    {len(pm)} records")
        pubmed_all.extend(pm)

        ep = fetch_epmc(query, args.max_per_query)
        print(f"         EuropePMC: {len(ep)} records")
        epmc_all.extend(ep)

    pubmed_deduped = deduplicate(pubmed_all)
    epmc_deduped = deduplicate(epmc_all)
    combined = deduplicate(pubmed_deduped + epmc_deduped)

    pubmed_out = OUTPUT_DIR / "pubmed_raw.json"
    epmc_out = OUTPUT_DIR / "epmc_raw.json"

    with open(pubmed_out, "w", encoding="utf-8") as f:
        json.dump(pubmed_deduped, f, indent=2, ensure_ascii=False)

    with open(epmc_out, "w", encoding="utf-8") as f:
        json.dump(epmc_deduped, f, indent=2, ensure_ascii=False)

    print("\n" + "=" * 60)
    print(f"PubMed unique:    {len(pubmed_deduped):>6,}")
    print(f"EuropePMC unique: {len(epmc_deduped):>6,}")
    print(f"Combined unique:  {len(combined):>6,}")
    print("=" * 60)

    if len(combined) < 1000:
        print(
            f"\nWARNING: {len(combined):,} papers collected — below 1000 target. "
            f"Increase --max-per-query or add more queries."
        )
    else:
        print(f"\nTarget met: {len(combined):,} unique papers collected.")

    print(f"\nSaved to:\n  {pubmed_out}\n  {epmc_out}")


if __name__ == "__main__":
    main()
