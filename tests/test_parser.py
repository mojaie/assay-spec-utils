import copy
from pathlib import Path
import tempfile
import unittest

from pydantic import ValidationError

from assay_spec_utils.model import AssayProtocol, AssayTemplates, AssayAttributes
from assay_spec_utils.parser import *
from assay_spec_utils.parser import _resolve_attributes

EXAMPLE_DIR = Path(__file__).resolve().parents[1] / "example"


def _attributes():
    attrs = [
        dict(attributeId="A1",
             terms=[("t:a1", "term a1")],
             parameters=[("concentration", "NaCl", 100, "mM")]),
        dict(attributeId="A2",
             terms=[("t:a2", "term a2")]),
    ]
    spec = AssayAttributes(items=attrs).model_dump()
    return {a["attributeId"]: a for a in spec["items"]}


def _templates():
    tmpls = [
        dict(
            templateId="TMPL",
            targets=[dict(accessionId="P00001")],
            attributes=["A1"],
            terms=[("t:tmpl", "template term")],
            parameters=[("time", "treatment", 1, "hour"),
                        ("volume", "total", 10, "uL")],
            readouts=[
                dict(readoutId="r1", terms=[("t:r1", "readout 1")]),
                dict(readoutId="r2", readoutMode="activation",
                     targets=[dict(accessionId="P00002")],
                     attributes=["A2"]),
            ]
        )
    ]
    spec = AssayTemplates(items=tmpls).model_dump()
    return {t["templateId"]: t for t in spec["items"]}


def _protocol(protocol_id="PROT", **fields):
    spec = dict(
        templateId="TMPL",
        terms=[("t:prot", "protocol term")],
        parameters=[("time", "treatment", 2, "hour")],
        assays=[
            dict(
                assayId="primary",
                valueType="percentage",
                terms=[("t:assay", "assay term")],
                parameters=[("volume", "total", 20, "uL")],
                datasources=[
                    dict(sourceType="CSV", sourcePath="a.csv",
                         sampleIdColumn="CID", valueColumns=["v1", None]),
                ]
            ),
            dict(
                assayId="screen",
                datasources=[
                    dict(sessionId="AS-0001", layerIndices=[3, 4],
                         subExperiments=[dict(subExpKey="plate",
                                              subExpValues=["x", "y"])]),
                ]
            ),
        ]
    )
    spec.update(fields)
    p = AssayProtocol(**spec).model_dump()
    p["protocolId"] = protocol_id
    return p


class TestParseSpec(unittest.TestCase):
    def test_parse_spec(self):
        spec = parse_spec({"templateId": "TMPL"}, AssayProtocol)
        self.assertEqual(spec["templateId"], "TMPL")
        self.assertEqual(spec["assays"], [])

    def test_parse_spec_invalid(self):
        with self.assertRaises(ValidationError):
            parse_spec({"unknownField": 1}, AssayProtocol)

    def test_parse_spec_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "p.yaml"
            path.write_text("templateId: TMPL\nterms:\n  - [t:1, term1]\n")
            spec = parse_spec_file(path, AssayProtocol)
        self.assertEqual(spec["terms"], [("t:1", "term1")])


class TestLoadFiles(unittest.TestCase):
    def test_load_skips_underscore_files(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "b.yaml").write_text("templateId: T\n")
            (d / "a.yaml").write_text("{}\n")
            (d / "_draft.yaml").write_text("templateId: T\n")
            (d / "note.txt").write_text("not a spec\n")
            protocols = load_protocols(d)
        self.assertEqual([p["protocolId"] for p in protocols], ["a", "b"])

    def test_load_example_protocols(self):
        protocols = load_protocols(EXAMPLE_DIR / "protocols")
        self.assertEqual(len(protocols), 5)
        self.assertIn("p53_Mdm2_HTRF", [p["protocolId"] for p in protocols])

    def test_load_example_templates(self):
        templates = load_templates(EXAMPLE_DIR / "templates")
        self.assertEqual(len(templates), 13)
        self.assertIn("DUAL_GLO", templates)

    def test_load_example_attributes(self):
        attributes = load_attributes(EXAMPLE_DIR / "attributes")
        self.assertEqual(len(attributes), 48)
        self.assertIn("PBS", attributes)


class TestResolveAttributes(unittest.TestCase):
    def test_resolve_attributes(self):
        spec = dict(
            attributes=["A1", "A2"],
            terms=[("t:own", "own term")],
            parameters=[("pH", "buffer", 7.4, None)]
        )
        _resolve_attributes(spec, _attributes())
        self.assertEqual(
            spec["terms"],
            ["t:own", "t:a1", "attr:A1", "t:a2", "attr:A2"])
        self.assertEqual(
            spec["parameters"],
            [("pH", "buffer", 7.4, None),
             ("concentration", "NaCl", 100, "mM")])

    def test_unknown_attribute(self):
        spec = dict(attributes=["UNKNOWN"], terms=[], parameters=[])
        with self.assertRaises(KeyError):
            _resolve_attributes(spec, _attributes())


class TestGenerateTermDict(unittest.TestCase):
    def test_collects_all_levels(self):
        termdict = generate_term_dict(
            [_protocol()], _templates(), _attributes())
        self.assertEqual(termdict, {
            "t:a1": "term a1", "t:a2": "term a2",
            "t:tmpl": "template term", "t:r1": "readout 1",
            "t:prot": "protocol term", "t:assay": "assay term"
        })


class TestGenerateProtocols(unittest.TestCase):
    def test_template_and_protocol_terms(self):
        specs = generate_protocols([_protocol()], _templates(), _attributes())
        self.assertEqual(specs, [{
            "protocolId": "PROT",
            "terms": ["t:tmpl", "t:a1", "attr:A1", "attr:TMPL", "t:prot"]
        }])

    def test_without_template(self):
        p = _protocol(templateId=None, attributes=["A2"])
        specs = generate_protocols([p], {}, _attributes())
        self.assertEqual(specs[0]["terms"], ["t:prot", "t:a2", "attr:A2"])


class TestGenerateAssays(unittest.TestCase):
    def setUp(self):
        self.assays = generate_assays(
            [_protocol()], _templates(), _attributes())
        self.by_key = {(a["assayId"], a["readoutId"]): a for a in self.assays}

    def test_readout_assay_combinations(self):
        # valueColumns[1] is null, so (primary, r2) is not generated
        self.assertEqual(
            sorted(self.by_key),
            [("primary", "r1"), ("screen", "r1"), ("screen", "r2")])

    def test_terms_unique_and_sorted(self):
        a = self.by_key[("primary", "r1")]
        self.assertEqual(a["terms"], sorted([
            "attr:A1", "attr:TMPL", "t:a1", "t:assay", "t:prot", "t:r1",
            "t:tmpl"
        ]))
        b = self.by_key[("screen", "r2")]
        self.assertIn("attr:A2", b["terms"])
        self.assertNotIn("t:r1", b["terms"])

    def test_parameters_overridden_by_lower_level(self):
        a = self.by_key[("primary", "r1")]
        # Keyed by (category, name); the last one wins and keeps its position
        self.assertEqual(a["parameters"], [
            ("concentration", "NaCl", 100, "mM"),
            ("time", "treatment", 2, "hour"),
            ("volume", "total", 20, "uL"),
        ])
        b = self.by_key[("screen", "r1")]
        self.assertIn(("volume", "total", 10, "uL"), b["parameters"])

    def test_targets_and_readout_fields(self):
        r1 = self.by_key[("screen", "r1")]
        r2 = self.by_key[("screen", "r2")]
        self.assertEqual([t["accessionId"] for t in r1["targets"]], ["P00001"])
        self.assertEqual([t["accessionId"] for t in r2["targets"]], ["P00002"])
        self.assertEqual(r1["readoutMode"], "inhibition")
        self.assertEqual(r2["readoutMode"], "activation")

    def test_protocol_targets_override_template(self):
        p = _protocol(targets=[dict(accessionId="P00009")])
        assays = generate_assays([p], _templates(), _attributes())
        by_key = {(a["assayId"], a["readoutId"]): a for a in assays}
        self.assertEqual(
            by_key[("screen", "r1")]["targets"][0]["accessionId"], "P00009")
        # Readout-level targets take precedence
        self.assertEqual(
            by_key[("screen", "r2")]["targets"][0]["accessionId"], "P00002")

    def test_datasources(self):
        self.assertEqual(
            self.by_key[("primary", "r1")]["datasources"][0]["valueColumn"],
            "v1")
        ds = self.by_key[("screen", "r2")]["datasources"][0]
        self.assertEqual(ds["sessionId"], "AS-0001")
        self.assertEqual(ds["layerIndex"], 4)
        self.assertEqual(
            ds["subExperiment"], [{"subExpKey": "plate", "subExpValue": "y"}])

    def test_value_type(self):
        self.assertEqual(
            self.by_key[("primary", "r1")]["valueType"], "percentage")
        self.assertEqual(self.by_key[("screen", "r1")]["valueType"], "AC50")

    def test_protocol_readouts_replace_template(self):
        p = _protocol(readouts=[dict(readoutId="own")])
        assays = generate_assays([p], _templates(), _attributes())
        self.assertEqual({a["readoutId"] for a in assays}, {"own"})

    def test_assays_without_readouts(self):
        p = _protocol(templateId=None)
        with self.assertRaises(ValueError):
            generate_assays([p], {}, _attributes())

    def test_example(self):
        protocols = load_protocols(EXAMPLE_DIR / "protocols")
        templates = load_templates(EXAMPLE_DIR / "templates")
        attributes = load_attributes(EXAMPLE_DIR / "attributes")
        assays = generate_assays(protocols, templates, attributes)
        self.assertEqual(len(assays), 9)
        keys = {(a["protocolId"], a["assayId"], a["readoutId"])
                for a in assays}
        self.assertIn(("p53_Mdm2_HTRF", "dose", "ratio"), keys)

    @unittest.expectedFailure
    def test_templates_reusable_after_generate_protocols(self):
        # Known bug: templates are resolved in place, so terms are broken
        # when the same templates dict is reused.
        templates = _templates()
        generate_protocols(
            [_protocol()], templates, _attributes())
        assays = generate_assays(
            [_protocol()], templates, _attributes())
        self.assertEqual(
            {(a["assayId"], a["readoutId"]): a["terms"] for a in assays},
            {(a["assayId"], a["readoutId"]): a["terms"]
             for a in self.assays})


if __name__ == "__main__":
    unittest.main()
