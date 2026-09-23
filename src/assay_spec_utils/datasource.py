
import logging
import pandas as pd
from pathlib import Path

import requests

__all__ = [
    "fetch_uniprot_json",
    "fetch_chebi_json",
    "chebi_name",
    "uniprot_chebi_terms",
    "uniprot_go_terms"
]
logger = logging.getLogger(__name__)

# External resources

# UniProt REST API
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/"
# EBI-OLS4 Web API
EBI_BASE_URL = "https://www.ebi.ac.uk/ols4/api/ontologies/"
# PUG REST API (PubChem)
PUG_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"


def fetch_uniprot_json(accession_id: str):
    r = requests.get(f"{UNIPROT_BASE_URL}{accession_id}.json")
    return r.json()


def fetch_chebi_json(obo_id: str):
    obo_num = obo_id.split(":")[1]
    query = f"{EBI_BASE_URL}chebi/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_{obo_num}"
    r = requests.get(query)
    return r.json()


def chebi_name(chebi_json: dict) -> str:
    """e.g. CHEBI:53438
    """
    return chebi_json["label"]


def uniprot_chebi_terms(uniprot_json: dict) -> tuple[dict, dict]:
    """Retrieve target GO and ChEBI terms from UniProt by accession ID
    """
    terms = []
    # target reactions
    for rcd in uniprot_json["comments"]:
        if rcd["commentType"] == "CATALYTIC ACTIVITY":
            if "reactionCrossReferences" not in rcd["reaction"]:
                continue
            for r in rcd["reaction"]["reactionCrossReferences"]:
                if r["database"] == "ChEBI":
                    terms.append(r["id"])
        elif rcd["commentType"] == "COFACTOR":
            for r in rcd["cofactors"]:
                cof = r["cofactorCrossReference"]
                if cof["database"] == "ChEBI":
                    terms.append(cof["id"])
    return terms


def uniprot_go_terms(uniprot_json: dict) -> tuple[dict, dict]:
    """Retrieve target GO and ChEBI terms from UniProt by accession ID
    """
    terms = {"Function": [], "Process": [], "Component": []}
    term_name = {}  # term ID => term name
    # target GO terms
    for rcd in uniprot_json["uniProtKBCrossReferences"]:
        if rcd["database"] != "GO":
            continue
        if rcd["properties"][0]["key"] != "GoTerm":
            raise ValueError("Unexpected format in uniProtKBCrossReferences")
        r = rcd["properties"][0]["value"]
        got, term = r.split(":")[:2]
        gotype = {"F": "Function", "P": "Process", "C": "Component"}[got]
        terms[gotype].append(rcd["id"])
        if rcd["id"] not in term_name:
            term_name[rcd["id"]] = " ".join(term.splitlines())
    return terms, term_name


def pubchem_assay(aid: str):
    query = f"{PUG_BASE_URL}/assay/aid/{aid}/concise/CSV"
    res = requests.get(query).json()
    return res["label"]


def load_table(spec, base_dir: Path, **pd_kwargs):
    """load dataset from a local CSV file"""
    assert spec["sourceType"] == "CSV"
    df = pd.read_csv(base_dir / spec["sourcePath"], **pd_kwargs)
    df.set_index(spec["sampleIdColumn"])
    # TODO: string to nan
    return df[spec["column"]]