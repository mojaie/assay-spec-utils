import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from assay_spec_utils.run import *
from tests.fixtures import CHEBI_LABELS, chebi_json, uniprot_json

TEMPLATES_YAML = """\
items:
  - templateId: TMPL
    targets:
      - accessionId: P00001
    readouts:
      - readoutId: r1
        targets:
          - accessionId: P00002
"""

PROTOCOL_YAML = """\
templateId: TMPL
targets:
  - accessionId: P00003
terms:
  - [t:prot, protocol term]
attributes:
  - A1
assays:
  - assayId: dose
    datasources:
      - sourceType: CSV
        sourcePath: a.csv
        sampleIdColumn: CID
        valueColumns: [IC50]
"""

ATTRIBUTES_YAML = """\
items:
  - attributeId: A1
    terms:
      - [t:a1, attribute term]
"""


class RunTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        base = Path(self._tmp.name)
        self.src = base / "src"
        self.dest = base / "dest"
        for d in ["protocols", "templates", "attributes"]:
            (self.src / d).mkdir(parents=True)
        for d in ["uniprot", "chebi", "processed"]:
            (self.dest / d).mkdir(parents=True)
        (self.src / "templates" / "t.yaml").write_text(TEMPLATES_YAML)
        (self.src / "protocols" / "PROT.yaml").write_text(PROTOCOL_YAML)
        (self.src / "attributes" / "a.yaml").write_text(ATTRIBUTES_YAML)

    def tearDown(self):
        self._tmp.cleanup()

    def _write_cache(self, accession_id):
        with open(self.dest / "uniprot" / f"{accession_id}.json", "w") as f:
            json.dump(uniprot_json(accession_id), f)
        for cid, label in CHEBI_LABELS.items():
            path = self.dest / "chebi" / f"{cid.split(':')[1]}.json"
            with open(path, "w") as f:
                json.dump(chebi_json(label), f)


class TestFetchTargetInfo(RunTestCase):
    def _fetch(self, **kwargs):
        with mock.patch(
                "assay_spec_utils.run.fetch_uniprot_json",
                side_effect=lambda aid: uniprot_json(
                    aid, reviewed=(aid != "P00002"))) as fu, \
             mock.patch(
                "assay_spec_utils.run.fetch_chebi_json",
                side_effect=lambda cid: chebi_json(CHEBI_LABELS[cid])) as fc:
            fetch_target_info(self.src, self.dest, **kwargs)
        return fu, fc

    def test_fetch(self):
        fu, fc = self._fetch()
        self.assertEqual(sorted(c.args[0] for c in fu.call_args_list),
                         ["P00001", "P00002", "P00003"])
        # Unreviewed entries are skipped
        self.assertEqual(
            sorted(p.name for p in (self.dest / "uniprot").iterdir()),
            ["P00001.json", "P00003.json"])
        self.assertEqual(
            sorted(p.name for p in (self.dest / "chebi").iterdir()),
            ["15377.json", "29105.json"])
        # ChEBI files already downloaded are not fetched again
        self.assertEqual(fc.call_count, 2)

    def test_cache(self):
        self._write_cache("P00001")
        fu, _ = self._fetch()
        self.assertNotIn("P00001", [c.args[0] for c in fu.call_args_list])

    def test_force_update(self):
        self._write_cache("P00001")
        fu, fc = self._fetch(force_update=True)
        self.assertIn("P00001", [c.args[0] for c in fu.call_args_list])
        self.assertEqual(fc.call_count, 4)


class TestTargetTerms(RunTestCase):
    def test_target_terms(self):
        self._write_cache("P00001")
        targets, termdict = target_terms(self.dest)
        self.assertEqual(targets, {"P00001": {
            "Function": ["GO:0003677"],
            "Process": ["GO:0006915"],
            "Component": ["GO:0005634"],
            "ChEBI": ["CHEBI:15377", "CHEBI:29105"]
        }})
        self.assertEqual(termdict["GO:0006915"], "apoptotic process")
        self.assertEqual(termdict["CHEBI:29105"], "zinc(2+)")

    def test_empty(self):
        self.assertEqual(target_terms(self.dest), ({}, {}))


class TestProcessAll(RunTestCase):
    def test_process_all(self):
        self._write_cache("P00003")
        process_all(self.src, self.dest)
        with gzip.open(self.dest / "processed" / "assays.json.gz", "rt") as f:
            result = json.load(f)
        self.assertEqual(set(result), {"meta", "assays", "targets", "terms"})
        self.assertIn("created", result["meta"])
        self.assertEqual(len(result["assays"]), 1)
        assay = result["assays"][0]
        self.assertEqual(
            (assay["protocolId"], assay["assayId"], assay["readoutId"]),
            ("PROT", "dose", "r1"))
        self.assertEqual(assay["targets"][0]["accessionId"], "P00002")
        self.assertEqual(list(result["targets"]), ["P00003"])
        self.assertEqual(result["terms"]["t:prot"], "protocol term")
        self.assertEqual(result["terms"]["t:a1"], "attribute term")
        self.assertEqual(result["terms"]["CHEBI:15377"], "water")


if __name__ == "__main__":
    unittest.main()
