import os
from pathlib import Path

from assay_spec_utils.parser import *
from assay_spec_utils.parser import _resolve_attributes
BASE_DIR = Path("./example")


def test_load_protocols():
    protocols = load_protocols(BASE_DIR / "protocols")
    assert len(protocols) == 2


def test_load_templates():
    templates = load_templates(BASE_DIR / "templates")
    assert len(templates) == 9


def test_load_attributes():
    attributes = load_attributes(BASE_DIR / "attributes")
    assert len(attributes) == 27


def test_resolve_attributes():
    attrs = dict(
        test1=dict(
            terms=["t1 hogehoge", "t2 fuga", "t3 piyo piyo"],
            parameters=dict(p1=10, p2="value", p3=1.23)
        )
    )
    spec = dict(
        attributes=["test1"],
        terms=["t4 hoge"],
        parameters=dict(p4="hoge")
    )
    _resolve_attributes(spec, attrs)
    assert spec["terms"] == ["t4", "t1", "t2", "t3"]
    assert spec["parameters"] == dict(p1=10, p2="value", p3=1.23, p4="hoge")
