
from datetime import datetime
import gzip
import json
import logging
from pathlib import Path

from assay_spec_utils.parser import *
from assay_spec_utils.datasource import *

__all__ = [
    "fetch_target_info",
    "target_terms",
    "process_all"
]
logger = logging.getLogger(__name__)


def _fetch_targets(
            spec: dict, dest_uniprot_dir: Path, dest_chebi_dir: Path,
            force_update=False
        ) -> None:
    for target in spec["targets"]:
        aid = target["accessionId"]
        dest_uniprot_file = dest_uniprot_dir / f"{aid}.json"
        if dest_uniprot_file.exists() and not force_update:
            continue
        uj = fetch_uniprot_json(aid)
        if uj["entryType"] != "UniProtKB reviewed (Swiss-Prot)":
            logger.warning(f"Skipped {aid}: invalid entryType {uj["entryType"]}")
            continue
        with open(dest_uniprot_file, "wt", encoding='UTF-8') as f:
            json.dump(uj, f, indent=2)
        cterms = uniprot_chebi_terms(uj)
        for ct in cterms:
            dest_chebi_file = dest_chebi_dir / f"{ct.split(":")[1]}.json"
            if dest_chebi_file.exists() and not force_update:
                continue
            cj = fetch_chebi_json(ct)
            with open(dest_chebi_file, "wt", encoding='UTF-8') as f:
                json.dump(cj, f, indent=2)


def fetch_target_info(
            src_dir: Path, dest_dir: Path,
            src_protocol_dir="protocols", src_templates_dir="templates",
            dest_uniprot_dir="uniprot", dest_chebi_dir="chebi",
            force_update=False
        ):
    protocols = load_protocols(src_dir / src_protocol_dir)
    templates = load_templates(src_dir / src_templates_dir)
    tdir = dest_dir / dest_uniprot_dir
    cdir = dest_dir / dest_chebi_dir
    # TODO: ncRNA, unknown gene
    for tid, tmpl in templates.items():
        logger.info(f"Targets: {tid}")
        _fetch_targets(tmpl, tdir, cdir, force_update)
        for readout in tmpl["readouts"]:
            _fetch_targets(readout, tdir, cdir, force_update)
    for protocol in protocols:
        logger.info(f"Targets: {protocol['protocolId']}")
        _fetch_targets(protocol, tdir, cdir, force_update)
        for readout in protocol["readouts"]:
            _fetch_targets(readout, tdir, cdir, force_update)
    logger.info("Done.")
    return


def target_terms(
            base_dir: Path,
            uniprot_dir="uniprot", chebi_dir="chebi",
        ):
    termdict = {}
    targets = {}
    for upath in sorted((base_dir / uniprot_dir).glob("*.json")):
        with open(upath, "rt", encoding='UTF-8') as f:
            uj = json.load(f)
        tterms, tgdict = uniprot_go_terms(uj)
        targets[upath.stem] = tterms
        cterms = uniprot_chebi_terms(uj)
        targets[upath.stem]["ChEBI"] = cterms
        cdict = {}
        for ct in cterms:
            cpath = base_dir / chebi_dir / f"{ct.split(":")[1]}.json"
            with open(cpath, "rt", encoding='UTF-8') as f:
                cj = json.load(f)
            cdict[ct] = chebi_name(cj)
        termdict.update(tgdict)
        termdict.update(cdict)
    return targets, termdict


def process_all(
            src_dir: Path, dest_dir: Path,
            src_protocol_dir="protocols", src_templates_dir="templates",
            src_attributes_dir="attributes", dest_processed_dir="processed",
            dest_assay_file="assays.json.gz"
        ):
    logger.info("Loading specification files...")
    protocols = load_protocols(src_dir / src_protocol_dir)
    templates = load_templates(src_dir / src_templates_dir)
    attributes = load_attributes(src_dir / src_attributes_dir)

    logger.info("Creating a term dictionary...")
    termdict = generate_term_dict(protocols, templates, attributes)

    logger.info("Generating target information...")
    targets, targetdict = target_terms(dest_dir)
    termdict.update(targetdict)

    logger.info("Generating assays...")
    assays = generate_assays(protocols, templates, attributes)
    assay_json = {
        "meta": {
            "created": datetime.now().isoformat(timespec="seconds")
        },
        "assays": assays,
        "targets": targets,
        "terms": termdict
    }
    assay_dest = dest_dir / dest_processed_dir / dest_assay_file
    with gzip.open(assay_dest, "wt", encoding='UTF-8') as f:
        json.dump(assay_json, f, indent=2)
    logger.info(f"Saved: {assay_dest}")
