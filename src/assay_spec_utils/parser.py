
import logging
from pathlib import Path
import pickle
import re
import yaml

from assay_spec_utils.datasource import uniprot_target_terms
from assay_spec_utils.model import AssayProtocol, AssayTemplates, AssayAttributes, AssaySpec

__all__ = [
    "parse_spec",
    "parse_spec_file",
    "load_protocols",
    "load_templates",
    "load_attributes",
    "fetch_target_terms",
    "generate_term_dict",
    "generate_assays"
]

logger = logging.getLogger(__name__)


def parse_spec(spec, validator, **kwargs):
    """Validate and parse dict-like assay spec data
    """
    return validator(**spec).model_dump(**kwargs)


def parse_spec_file(path, validator, **kwargs):
    """Validate and parse assay spec file"""
    with open(path) as f:
        spec = yaml.safe_load(f)
    return parse_spec(spec, validator, **kwargs)


def load_protocols(protocols_dir: Path, **kwargs) -> list:
    """load assay protocol files"""
    protocols = []
    for fpath in sorted(protocols_dir.glob("*.yaml")):
        if fpath.name.startswith("_"):
            continue  # Example or draft
        logger.info(f"Protocol: {fpath.name}")
        spec = parse_spec_file(fpath, AssayProtocol, **kwargs)
        spec["protocolId"] = fpath.stem
        protocols.append(spec)
    logger.info(f"Done.")
    return protocols


def load_templates(templates_dir: Path, **kwargs) -> dict:
    """load protocol template files"""
    templates = {}
    for fpath in sorted(templates_dir.glob("*.yaml")):
        if fpath.name.startswith("_"):
            continue  # Example or draft
        logger.info(f"Templates: {fpath.name}")
        for tmpl in parse_spec_file(fpath, AssayTemplates, **kwargs)["items"]:
            templates[tmpl["templateId"]] = tmpl
    logger.info(f"Done.")
    return templates


def load_attributes(attrs_dir: Path, **kwargs) -> dict:
    """load protocol template files"""
    attributes = {}
    for fpath in sorted(attrs_dir.glob("*.yaml")):
        if fpath.name.startswith("_"):
            continue  # Example or draft
        logger.info(f"Attributes: {fpath.name}")
        for attr in parse_spec_file(fpath, AssayAttributes, **kwargs)["items"]:
            attributes[attr["attributeId"]] = attr
    logger.info(f"Done.")
    return attributes


def _fetch_targets(spec, target_term, term_dict) -> None:
    for target in spec["targets"]:
        tgtm, tmnm = uniprot_target_terms(target["accessionId"])
        target_term[target["accessionId"]] = tgtm
        term_dict.update(tmnm)


def fetch_target_terms(protocols: list, templates: dict) -> tuple[dict, dict]:
    """Fetch target terms from UniProt"""
    # TODO: ncRNA, unknown gene
    target_term = {}  # UniProtID => {GOtype => [GOterms]}
    term_dict = {}  # GOterm => GOname
    for tid, tmpl in templates.items():
        logger.info(f"Template: {tid}")
        _fetch_targets(tmpl, target_term, term_dict)
        for readout in tmpl["readouts"]:
            _fetch_targets(readout, target_term, term_dict)
    for protocol in protocols:
        logger.info(f"Protocol: {protocol['protocolId']}")
        _fetch_targets(protocol, target_term, term_dict)
        for readout in protocol["readouts"]:
            _fetch_targets(readout, target_term, term_dict)
    logger.info("Done.")
    return target_term, term_dict


def _update_term_dict(spec, term_dict) -> None:
    for term, name in spec["terms"]:
        term_dict[term] = name


def generate_term_dict(
        protocols: list, templates: dict, attributes: dict) -> dict:
    """Generate dict of terms that appear in all specification files"""
    term_dict = {}
    for _, attr in attributes.items():
        _update_term_dict(attr, term_dict)
    for _, tmpl in templates.items():
        _update_term_dict(tmpl, term_dict)
        for readout in tmpl["readouts"]:
            _update_term_dict(readout, term_dict)
    for protocol in protocols:
        _update_term_dict(protocol, term_dict)
        for readout in protocol["readouts"]:
            _update_term_dict(readout, term_dict)
        for assay in protocol["assays"]:
            _update_term_dict(assay, term_dict)
    return term_dict


def _resolve_attributes(spec, gattrs) -> None:
    spec["terms"] = [t[0] for t in spec["terms"]]
    for attrid in spec["attributes"]:
        attr = gattrs[attrid]
        # Should be unique and sorted later
        atms = [t[0] for t in attr["terms"]]
        spec["terms"] = [*spec["terms"], *atms]
        spec["parameters"] = [*spec["parameters"], *attr["parameters"]]


def generate_assays(
        protocols: list, templates: dict,
        attributes: dict, target_term: dict) -> list:
    """Generate assay metadata for each datasets.
    """
    rcds = []
    templates_resolved = {}
    templates_readouts = {}
    for protocol in protocols:
        logger.info(f"Protocol: {protocol['protocolId']}")

        # Apply template
        tmpl_id = protocol["templateId"]
        if tmpl_id is None:
            spec = {
                "targets": [],
                "terms": [],
                "parameters": []
            }
            readouts = []
        else:
            if tmpl_id not in templates_resolved:
                tmpl = templates[tmpl_id]
                _resolve_attributes(tmpl, attributes)
                templates_resolved[tmpl_id] = {
                    "targets": tmpl["targets"],
                    "terms": tmpl["terms"],
                    "parameters": tmpl["parameters"]
                }
                templates_readouts[tmpl_id] = tmpl["readouts"]
            spec = pickle.loads(pickle.dumps(templates_resolved[tmpl_id]))
            readouts = pickle.loads(pickle.dumps(templates_readouts[tmpl_id]))

        # Override by protocol-level fields
        if protocol["targets"]:
            spec["targets"] = protocol["targets"]
        _resolve_attributes(protocol, attributes)
        spec["terms"].extend(protocol["terms"])
        spec["parameters"].extend(protocol["parameters"])
        spec["protocolId"] = protocol["protocolId"]
        if protocol["readouts"]:
            readouts = protocol["readouts"]
        if not readouts and protocol["assays"]:
            raise ValueError("Assays defined but no readouts found.")

        # Readouts
        for ridx, readout in enumerate(readouts):
            sp = pickle.loads(pickle.dumps(spec))
            ro = pickle.loads(pickle.dumps(readout))
            # Resolve targets
            if ro["targets"]:
                sp["targets"] = ro["targets"]
            for tg in sp["targets"]:
                tt = target_term[tg["accessionId"]]
                sp["terms"].extend(tt["Function"])
                sp["terms"].extend(tt["Process"])
                sp["terms"].extend(tt["Component"])
                sp["terms"].extend(tt["ChEBI"])

            # Override by readout-level fields
            _resolve_attributes(ro, attributes)
            sp["terms"].extend(ro["terms"])
            sp["parameters"].extend(ro["parameters"])
            sp["readoutId"] = ro["readoutId"]
            sp["readoutType"] = ro["readoutType"]

            # Assays
            for assay in protocol["assays"]:
                #Data
                data = []
                for d in assay["datasources"]:
                    newd = {"sourceType": d["sourceType"]}
                    if d["sourceType"] == "Screener":
                        if d["layerIndices"][ridx] is None:
                            continue
                        newd["sessionId"] = d["sessionId"]
                        newd["layerIndex"] = d["layerIndices"][ridx]
                        if d["subExperiments"]:
                            newd["subExperiment"] = []
                            for i, subex in enumerate(d["subExperiments"]):
                                newd["subExperiment"].append({
                                    "subExpKey": subex["subExpKey"],
                                    "subExpValue": subex["subExpValues"][ridx],
                                })
                    elif d["sourceType"] == "CSV":
                        if d["valueColumns"][ridx] is None:
                            continue
                        newd["sourcePath"] = d["sourcePath"]
                        newd["sampleIdColumn"] = d["sampleIdColumn"]
                        newd["valueColumn"] = d["valueColumns"][ridx]
                    newd["sampleMapping"] = d["sampleMapping"]
                    data.append(newd)
                if not data:
                    continue
                s = pickle.loads(pickle.dumps(sp))
                asy = pickle.loads(pickle.dumps(assay))
                s["datasources"] = data

                # Override by assay-level fields
                _resolve_attributes(asy, attributes)
                s["terms"].extend(asy["terms"])
                s["parameters"].extend(asy["parameters"])
                # Make terms unique and sorted
                s["terms"] = sorted(set(s["terms"]))
                # Override parameters
                newps = []
                newpks = set()
                for p in reversed(s["parameters"]):
                    if (p[0], p[1]) not in newpks:
                        newpks.add((p[0], p[1]))
                        newps.append(p)
                s["parameters"] = list(reversed(newps))
                s["assayId"] = asy["assayId"]
                s["valueType"] = asy["valueType"]

                # Validation
                s = parse_spec(s, AssaySpec)
                rcds.append(s)
    logger.info("Done.")
    return rcds
