
import pandas as pd
from pathlib import Path

import requests


# External resources

# UniProt REST API
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/"
# EBI-OLS4 Web API
EBI_BASE_URL = "https://www.ebi.ac.uk/ols4/api/ontologies/"
# PUG REST API (PubChem)
PUG_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"


def uniprot_target_terms(accession_id: str) -> tuple[dict, dict]:
    """Retrieve target GO and ChEBI terms from UniProt by accession ID
    """
    r = requests.get(f"{UNIPROT_BASE_URL}{accession_id}.json")
    res = r.json()

    termids = {"Function": [], "Process": [], "Component": [], "ChEBI": []}
    term_name = {}  # term ID => term name
    # target reactions
    for rcd in res["comments"]:
        if rcd["commentType"] == "CATALYTIC ACTIVITY":
            if "reactionCrossReferences" not in rcd["reaction"]:
                continue
            for r in rcd["reaction"]["reactionCrossReferences"]:
                if r["database"] == "ChEBI":
                    termids["ChEBI"].append(r["id"])
                    if r["id"] not in term_name:
                        term_name[r["id"]] = chebi_name(r["id"])
        elif rcd["commentType"] == "COFACTOR":
            for r in rcd["cofactors"]:
                cof = r["cofactorCrossReference"]
                if cof["database"] == "ChEBI":
                    termids["ChEBI"].append(cof["id"])
                    if cof["id"] not in term_name:
                        term_name[cof["id"]] = chebi_name(cof["id"])
    # target GO terms
    for rcd in res["uniProtKBCrossReferences"]:
        if rcd["database"] != "GO":
            continue
        if rcd["properties"][0]["key"] != "GoTerm":
            raise ValueError("Unexpected format in uniProtKBCrossReferences")
        r = rcd["properties"][0]["value"]
        got, term = r.split(":")[:2]
        gotype = {"F": "Function", "P": "Process", "C": "Component"}[got]
        termids[gotype].append(rcd["id"])
        if rcd["id"] not in term_name:
            term_name[rcd["id"]] = " ".join(term.splitlines())
    return termids, term_name


def chebi_name(obo_id: str) -> str:
    """find chebi name label by obo_id (e.g. CHEBI:53438)
    """
    obo_num = obo_id.split(":")[1]
    query = f"{EBI_BASE_URL}chebi/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCHEBI_{obo_num}"
    res = requests.get(query).json()
    return res["label"]


def pubchem_assay(aid: str):
    query = f"{PUB_BASE_URL}/assay/aid/{aid}/concise/CSV"
    res = requests.get(query).json()
    return res["label"]


def load_table(spec, base_dir: Path, **pd_kwargs):
    """load dataset from a local CSV file"""
    assert spec["sourceType"] == "CSV"
    df = pd.read_csv(base_dir / spec["sourcePath"], **pd_kwargs)
    df.set_index(spec["sampleIdColumn"])
    # TODO: string to nan
    return df[spec["column"]]