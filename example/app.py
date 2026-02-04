
# Example of assay metadata input form (Streamlit web application)

from io import StringIO
from pathlib import Path
import re
import uuid

import numpy as np
import pandas as pd
import streamlit as st
import yaml

from assay_spec_utils.model import AssayTemplates, AssayAttributes, AssayProtocol
from assay_spec_utils.parser import parse_spec_file
from assay_spec_utils.util import is_convertible_to_int, is_convertible_to_float


# Variables

TEMPLATES_DIR = Path("./templates")
ATTRIBUTES_DIR = Path("./attributes")

# yaml literal format
# https://stackoverflow.com/questions/6432605/any-yaml-libraries-in-python-that-support-dumping-of-long-strings-as-block-liter

class literal_str(str): pass

def literal_str_representer(dumper, data):
    return dumper.represent_scalar(u'tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal_str, literal_str_representer)


# TODO: min and max values
PARAM_TYPE_OPTIONS = [
    dict(name="temperature", type=int, units=["degC"]),
    dict(name="time", type=int, units=["sec", "min", "hour", "day"], min_value=0),
    dict(
        name="concentration",
        type=float,
        units=[
            "M", "mM", "uM", "nM", "pM",
            "g/L", "mg/L", "ug/L", "ng/L", "v/v%", "ratio", "dilution"
        ], min_value=0),
    dict(name="volume", type=float, units=["L", "mL", "uL", "nL", "pL"], min_value=0),
    dict(name="pH", type=float, units=None, min_value=0, max_value=14, step=0.1),
    dict(name="wavelength", type=int, units=["nm"], min_value=200, max_value=900),
    dict(name="cellcount", type=int, units=["cells/well"], min_value=0, step=100),
    dict(name="sequence", type=str, units=None),
    dict(name="position", type=str, units=None),
    dict(name="mutation", type=str, units=None),
    dict(name="integer", type=int, units=None),
    dict(name="float", type=float, units=None),
    dict(name="string", type=str, units=None)
]
PARAM_TYPE_DICT = {r["name"]: r for r in PARAM_TYPE_OPTIONS}


# Main logics

def template_options() -> list:
    templates = []
    files = ["binding", "enzyme", "reporter", "viability"]
    for f in files:
        templates.extend(parse_spec_file(
            TEMPLATES_DIR / f"{f}.yaml", AssayTemplates, exclude_unset=True)["items"])
    return templates


def attribute_options(attr: str) -> dict:
    return parse_spec_file(
        ATTRIBUTES_DIR / f"{attr}.yaml", AssayAttributes)["items"]


def load_template(saved_file, tmpl):
    if saved_file is not None:
        st.session_state["data"] = yaml.safe_load(saved_file)
    elif tmpl is not None:
        st.session_state["data"] = {
            "assayProtocolVersion": "1.0",
            "meta": {
                "description": ""
            },
            "targets": [],
            "attributes": tmpl.get("attributes", []),
            "terms": tmpl.get("terms", []),
            "parameters": tmpl.get("parameters", []),
            "readouts": tmpl["readouts"],
            "assays": []
        }
    else:
        return


def add_target():
    st.session_state["data"]["targets"].append(
        {"sourceType": "UniProt", "accessionId": ""})
    clear_generated()


def update_targets():
    new_targets = []
    for i, rcd in enumerate(st.session_state["data"]["targets"]):
        tid = st.session_state[f"tid_{i}"]
        new_targets.append({"sourceType": "UniProt", "accessionId": tid})
    st.session_state["data"]["targets"] = new_targets
    clear_generated()


def remove_target(idx):
    st.session_state["data"]["targets"].pop(idx)
    clear_generated()


def update_attributes():
    st.session_state["data"]["attributes"] = [
        st.session_state["cat"], st.session_state["plate"],
        *st.session_state["atag"], *st.session_state["reagent"]]
    clear_generated()


def add_parameter():
    st.session_state["data"]["parameters"].append(
        ["string", "New parameter", "", None])
    clear_generated()


def update_parameters():
    new_params = []
    for i, rcd in enumerate(st.session_state["data"]["parameters"]):
        pcat = st.session_state[f"pcat_{i}"]
        p = PARAM_TYPE_DICT[pcat]
        pname = st.session_state[f"pname_{i}"]
        pval = st.session_state[f"pval_{i}"]
        punit = st.session_state.get(f"punit_{i}")
        # set default values and units for the selected category
        if p["type"] is int and not is_convertible_to_int(pval):
            pval = 0
        elif p["type"] is float and not is_convertible_to_float(pval):
            pval = 0.0
        if p["units"] is None:
            punit = None
        elif punit not in p["units"]:
            punit = p["units"][0]
        new_params.append([pcat, pname, pval, punit])
    st.session_state["data"]["parameters"] = new_params
    clear_generated()


def remove_parameter(idx):
    st.session_state["data"]["parameters"].pop(idx)
    clear_generated()


def add_assay():
    st.session_state["data"]["assays"].append({
        "assayId": "primary",
        "valueType": "percentage",
        "datasources": [{
            "sourceType": "Screener",
            "sessionId": "",
            "layerIndices": [None for _ in st.session_state["data"]["readouts"]]
        }]
    })
    clear_generated()


def update_assays():
    new_assays = []
    for i, rcd in enumerate(st.session_state["data"]["assays"]):
        apid = st.session_state[f"apid_{i}"]
        aid = st.session_state[f"aid_{i}"]
        avtype = st.session_state[f"avtype_{i}"]
        asids = st.session_state[f"asid_{i}"]
        pval = st.session_state[f"pval_{i}"]
        assay = {
            "assayId": aid, 
            "valueType": avtype,
            "datasources": []
        }
        alys = []
        for j, readout in enumerate(st.session_state["data"]["readouts"]):
            alys.append(st.session_state[f"aly_{i}_{j}"])
        for sid in asids.split(","):
            assay["datasources"].append({
                "sourceType": "Screener",
                "sessionId": sid.strip(),
                "layerIndices": alys
            })
        new_assays.append(assay)
    st.session_state["data"]["assays"] = new_assays
    clear_generated()


def remove_assay(idx):
    st.session_state["data"]["assays"].pop(idx)
    clear_generated()


def update_description():
    st.session_state["data"]["meta"]["description"] = literal_str(st.session_state["description"])
    clear_generated()


def customize_yaml_style(yaml_dump, layer_cnt):
    # YAML style hack
    idt4list = "    - (.+?)"
    idt2list = "  - (.+?)"
    idt2list4 = "\n".join(idt2list for _ in range(4))
    idt4list4 = "\n".join(idt4list for _ in range(4))
    idt2list2 = "\n".join(idt2list for _ in range(2))
    idt4list2 = "\n".join(idt4list for _ in range(2))
    idt4listn = "\n".join(idt4list for _ in range(layer_cnt))
    rep4 = ", ".join(fr"\{i}" for i in range(1, 5))
    rep2 = ", ".join(fr"\{i}" for i in range(1, 3))
    repn = ", ".join(fr"\{i}" for i in range(1, layer_cnt + 1))
    yaml_dump = re.sub(
        fr"- !!python/tuple\n{idt2list4}\n", fr"- [{rep4}]\n", yaml_dump)
    yaml_dump = re.sub(
        fr"- !!python/tuple\n{idt4list4}\n", fr"- [{rep4}]\n", yaml_dump)
    yaml_dump = re.sub(
        fr"- !!python/tuple\n{idt2list2}\n", fr"- [{rep2}]\n", yaml_dump)
    yaml_dump = re.sub(
        fr"- !!python/tuple\n{idt4list2}\n", fr"- [{rep2}]\n", yaml_dump)
    yaml_dump = re.sub(
        fr"  layerIndices:\n{idt4listn}\n", fr"  layerIndices: [{repn}]\n", yaml_dump)
    return yaml_dump


def clear_generated():
    st.session_state.generated = False


def set_generated():
    st.session_state.generated = True


# UI

def run():
    st.title("Assay metadata form example")

    st.markdown("""
- Enter your assay details, and download a machine-readable assay metadata in YAML format.
- **Do not expose this application to any untrusted network, including the Internet.**""")

    saved_file = st.file_uploader(
        "Choose a saved assay file (.yml, .yaml)", type=["yaml", "yml"])
    tmpl_ops = template_options()
    tmpl_map = {rcd["templateId"]: rcd for rcd in tmpl_ops}
    tmpl_id = st.selectbox(
        "or choose an assay template", tmpl_map.keys(), key="templates", 
        index=None, format_func=lambda k: tmpl_map[k]["name"])
    col1, col2 = st.columns([1, 3])
    with col1:
        st.button("Apply template", type="primary",
            on_click=load_template, args=[saved_file, tmpl_map.get(tmpl_id)])
    with col2:
        if st.session_state.get("data"):
            st.write("to discard changes and start from a loaded file or a new template")
    if not st.session_state.get("data"):
        return


    st.markdown("""## Assay design""")

    # Targets
    for i, rcd in enumerate(st.session_state["data"]["targets"]):
        col1, col2, col3 = st.columns([1, 1, 2])
        with col1:
            st.text_input(
                "Target UniProt Accession ID", value=rcd["accessionId"],
                key=f"tid_{i}", on_change=update_targets)
        with col2:
            st.button("🗑️", key=f"tdel_{i}", on_click=remove_target, args=[i])
    st.button("Add a target", on_click=add_target)

    # Attributes
    attr = st.session_state["data"]["attributes"]
    cat_ops = attribute_options("assay_category")
    cat_map = {rcd["attributeId"]: rcd for rcd in cat_ops}
    cat = list(set(attr) & set(cat_map.keys()))
    cat_idx = list(cat_map.keys()).index(cat[0]) if cat else 0
    st.selectbox("Assay category", cat_map.keys(),
        key="cat", index=cat_idx, format_func=lambda k: cat_map[k]["name"],
        on_change=update_attributes)
    plate_ops = attribute_options("microplate")
    plate_map = {rcd["attributeId"]: rcd for rcd in plate_ops}
    plate = list(set(attr) & set(plate_map.keys()))
    plate_idx = list(plate_map.keys()).index(plate[0]) if plate else 0
    st.selectbox("Microplate", plate_map.keys(),
        key="plate", index=plate_idx, format_func=lambda k: plate_map[k]["name"],
        on_change=update_attributes)
    atag_ops = attribute_options("affinity_tag")
    atag_map = {rcd["attributeId"]: rcd for rcd in atag_ops}
    atags = list(set(attr) & set(atag_map.keys()))
    st.multiselect("Affinity tags", atag_map.keys(),
        default=atags, key="atag", format_func=lambda k: atag_map[k]["name"],
        on_change=update_attributes)
    reagent_ops = attribute_options("reagent")
    reagent_map = {rcd["attributeId"]: rcd for rcd in reagent_ops}
    reagents = list(set(attr)& set(reagent_map.keys()))
    st.multiselect("Mediums/premixed reagents", reagent_map.keys(),
        default=reagents, key="reagent",
        format_func=lambda k: reagent_map[k]["name"],
        on_change=update_attributes)


    st.markdown("""## Assay conditions""")

    for i, rcd in enumerate(st.session_state["data"]["parameters"]):
        col1, col2, col3, col4, col5 = st.columns([3, 3, 3, 3, 1])
        ptype_opts = list(PARAM_TYPE_DICT.keys())
        vtype = PARAM_TYPE_DICT[rcd[0]]["type"]
        with col1:
            st.selectbox(
                "Type", ptype_opts, index=ptype_opts.index(rcd[0]),
                key=f"pcat_{i}", on_change=update_parameters)
        with col2:
            st.text_input(
                "Name", value=rcd[1], key=f"pname_{i}",
                on_change=update_parameters)
        with col3:
            if vtype is int:
                st.number_input(
                    "Value", value=int(rcd[2]), key=f"pval_{i}", step=1,
                    on_change=update_parameters)
            elif vtype is float:
                st.number_input(
                    "Value", value=float(rcd[2]), key=f"pval_{i}",
                    on_change=update_parameters)
            elif vtype is str:
                st.text_input(
                    "Value", value=str(rcd[2]), key=f"pval_{i}",
                    on_change=update_parameters)
        with col4:
            punit_opts = PARAM_TYPE_DICT[rcd[0]]["units"]
            if punit_opts is not None:
                st.selectbox(
                    "Unit", punit_opts, index=punit_opts.index(rcd[3]),
                    key=f"punit_{i}")
        with col5:
            st.button(
                "🗑️", key=f"pdel_{i}", on_click=remove_parameter, args=[i])
    st.button("Add a parameter", on_click=add_parameter)


    st.markdown("""## Protocol descriptions and comments""")

    default_pid = saved_file.name if saved_file is not None else "new_assay"
    protocol_id = st.text_input(
        "Assay protocol ID", value=default_pid, key="protocol_id")
    st.text_area(
            "Description", key="description",
            value=st.session_state["data"]["meta"]["description"],
            height=150, on_change=update_description,
            placeholder="""e.g.
- Protein sample details (affinitiy tags, mutations, isoform, sequence)
- Sample lot information
- Experimental condition
""")


    st.markdown("""## Data source""")

    vtype_opts = ["percentage", "AC50", "RZ-score"]
    for i, rcd in enumerate(st.session_state["data"]["assays"]):
        col1, col2 = st.columns(2)
        default_sessions = ",".join(d["sessionId"] for d in rcd["datasources"])
        # For simplicity, all sessions are assumed to have the same layer settings
        default_layers = rcd["datasources"][0]["layerIndices"]
        with col1:
            st.text_input(
                "Assay ID", value=f"{protocol_id}_", disabled=True,
                label_visibility="hidden", key=f"apid_{i}",
                on_change=update_assays)
        with col2:
            st.text_input(
                "Assay ID", value=rcd["assayId"], key=f"aid_{i}",
                on_change=update_assays)
        st.selectbox(
            "Value type", vtype_opts,
            index=vtype_opts.index(rcd["valueType"]), key=f"avtype_{i}",
            on_change=update_assays)
        st.text_input(
            "Analyzer Session IDs (Comma separated)",
            value=default_sessions, key=f"asid_{i}",
            on_change=update_assays)
        for j, readout in enumerate(st.session_state["data"]["readouts"]):
            st.number_input(
                f"Layer index of {readout["readoutId"]} (1-based)",
                min_value=1, max_value=20, step=1,
                key=f"aly_{i}_{j}", value=default_layers[j],
                on_change=update_assays)
    st.button("Add an assay", on_click=add_assay)


    st.markdown("""## Generate and export assay metadata""")

    st.button("Generate assay metadata", type="primary", on_click=set_generated)
    if not st.session_state.get("generated"):
        return
    data = AssayProtocol(**st.session_state["data"]).model_dump(exclude_unset=True)
    st.write(data)  # display data and check
    dp = customize_yaml_style(
        yaml.dump(data, sort_keys=False), len(st.session_state["data"]["readouts"]))
    st.download_button(
        "Export YAML file", dp, file_name=f"{protocol_id}.yaml",
        mime="application/yaml")


run()