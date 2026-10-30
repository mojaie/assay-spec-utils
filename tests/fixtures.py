"""Minimal API responses for tests (no network access)"""


def uniprot_json(accession_id="P00001", reviewed=True):
    entry_type = ("UniProtKB reviewed (Swiss-Prot)" if reviewed
                  else "UniProtKB unreviewed (TrEMBL)")
    return {
        "entryType": entry_type,
        "primaryAccession": accession_id,
        "comments": [
            {
                "commentType": "CATALYTIC ACTIVITY",
                "reaction": {
                    "reactionCrossReferences": [
                        {"database": "Rhea", "id": "RHEA:10000"},
                        {"database": "ChEBI", "id": "CHEBI:15377"},
                    ]
                }
            },
            {
                "commentType": "CATALYTIC ACTIVITY",
                "reaction": {"name": "no cross references"}
            },
            {
                "commentType": "COFACTOR",
                "cofactors": [
                    {"cofactorCrossReference": {
                        "database": "ChEBI", "id": "CHEBI:29105"}}
                ]
            },
            {"commentType": "FUNCTION"},
        ],
        "uniProtKBCrossReferences": [
            {"database": "GO", "id": "GO:0003677", "properties": [
                {"key": "GoTerm", "value": "F:DNA binding"}]},
            {"database": "GO", "id": "GO:0006915", "properties": [
                {"key": "GoTerm", "value": "P:apoptotic\nprocess"}]},
            {"database": "GO", "id": "GO:0005634", "properties": [
                {"key": "GoTerm", "value": "C:nucleus"}]},
            {"database": "PDB", "id": "1A1U", "properties": []},
        ]
    }


def chebi_json(label):
    return {"label": label}


CHEBI_LABELS = {"CHEBI:15377": "water", "CHEBI:29105": "zinc(2+)"}
