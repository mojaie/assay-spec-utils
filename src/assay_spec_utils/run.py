
from datetime import datetime
import gzip
import json
import logging
from pathlib import Path

from assayspec.parser import *

logger = logging.getLogger(__name__)

PROTOCOLS_DIR = "protocols"
TEMPLATES_DIR = "templates"
ATTRIBUTES_DIR = "attributes"
TARGET_FILE = "targets.json.gz"
ASSAY_FILE = "assays.json.gz"
TERMDICT_FILE = "termdict.json.gz"


def process_all(spec_base_dir: Path, output_dir: Path):
    logger.info("Loading specification files...")
    protocols = load_protocols(spec_base_dir / PROTOCOLS_DIR)
    templates = load_templates(spec_base_dir / TEMPLATES_DIR)
    attributes = load_attributes(spec_base_dir / ATTRIBUTES_DIR)

    target_dest = output_dir / TARGET_FILE
    if target_dest.exists():
        logger.info("Loading target file...")
        with gzip.open(target_dest, "rt", encoding='UTF-8') as f:
            target_json = json.load(f)
        logger.info(f"Loaded: {target_dest}")
        targets = target_json["targets"]
        target_dict = target_json["term_dict"]
    else:
        logger.info("Target file not found. Fetching targets...")
        targets, target_dict = fetch_target_terms(protocols, templates)
        target_json = {
            "meta": {
                "created": datetime.now().isoformat(timespec="seconds")
            },
            "targets": targets,
            "term_dict": target_dict
        }
        with gzip.open(target_dest, "wt", encoding='UTF-8') as f:
            json.dump(target_json, f, indent=2)
        logger.info(f"Saved: {target_dest}")

    logger.info("Creating a term dictionary...")
    term_dict = generate_term_dict(protocols, templates, attributes)
    term_dict.update(target_dict)

    logger.info("Generating assays...")
    assays = generate_assays(protocols, templates, attributes, targets)
    with gzip.open(output_dir / ASSAY_FILE, "wt", encoding='UTF-8') as f:
        json.dump(assays, f, indent=2)
    logger.info(f"Saved: {output_dir / ASSAY_FILE}")
    with gzip.open(output_dir / TERMDICT_FILE, "wt", encoding='UTF-8') as f:
        json.dump(term_dict, f, indent=2)
    logger.info(f"Saved: {output_dir / TERMDICT_FILE}")