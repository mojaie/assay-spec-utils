
from datetime import datetime
import gzip
import json
import logging
from pathlib import Path

from assay_spec_utils.parser import *

logger = logging.getLogger(__name__)

PROTOCOLS_DIR = "protocols"
TEMPLATES_DIR = "templates"
ATTRIBUTES_DIR = "attributes"
TARGET_FILE = "targets.json.gz"
ASSAY_FILE = "assays.json.gz"


def process_all(src_dir: Path, dest_dir: Path):
    logger.info("Loading specification files...")
    protocols = load_protocols(src_dir / PROTOCOLS_DIR)
    templates = load_templates(src_dir / TEMPLATES_DIR)
    attributes = load_attributes(src_dir / ATTRIBUTES_DIR)

    target_dest = dest_dir / TARGET_FILE
    if target_dest.exists():
        logger.info("Loading target file...")
        with gzip.open(target_dest, "rt", encoding='UTF-8') as f:
            target_json = json.load(f)
        logger.info(f"Loaded: {target_dest}")
        targets = target_json["targets"]
        target_terms = target_json["terms"]
    else:
        logger.info("Target file not found. Fetching targets...")
        targets, target_terms = fetch_target_terms(protocols, templates)
        target_json = {
            "meta": {
                "created": datetime.now().isoformat(timespec="seconds")
            },
            "targets": targets,
            "terms": target_terms
        }
        with gzip.open(target_dest, "wt", encoding='UTF-8') as f:
            json.dump(target_json, f, indent=2)
        logger.info(f"Saved: {target_dest}")

    logger.info("Creating a term dictionary...")
    terms = generate_term_dict(protocols, templates, attributes)
    terms.update(target_terms)

    logger.info("Generating assays...")
    assays = generate_assays(protocols, templates, attributes, targets)
    assay_json = {
        "meta": {
            "created": datetime.now().isoformat(timespec="seconds")
        },
        "assays": assays,
        "targets": targets,
        "terms": terms
    }
    with gzip.open(dest_dir / ASSAY_FILE, "wt", encoding='UTF-8') as f:
        json.dump(assay_json, f, indent=2)
    logger.info(f"Saved: {dest_dir / ASSAY_FILE}")
