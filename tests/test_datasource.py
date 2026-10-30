from pathlib import Path
import tempfile
import unittest
from unittest import mock

import httpx2

from assay_spec_utils import datasource
from assay_spec_utils.datasource import *
from tests.fixtures import uniprot_json


def _response(url, status_code=200, **kwargs):
    return httpx2.Response(
        status_code, request=httpx2.Request("GET", url), **kwargs)


class TestFetch(unittest.TestCase):
    def test_fetch_uniprot_json(self):
        url = f"{datasource.UNIPROT_BASE_URL}P04637.json"
        with mock.patch("httpx2.get",
                        return_value=_response(url, json={"a": 1})) as get:
            self.assertEqual(fetch_uniprot_json("P04637"), {"a": 1})
        get.assert_called_once_with(
            url, follow_redirects=True, timeout=datasource.HTTP_TIMEOUT)

    def test_fetch_chebi_json(self):
        with mock.patch("httpx2.get", side_effect=lambda url, **_:
                        _response(url, json={"label": "water"})) as get:
            self.assertEqual(fetch_chebi_json("CHEBI:15377"),
                             {"label": "water"})
        url = get.call_args.args[0]
        self.assertTrue(url.startswith(f"{datasource.EBI_BASE_URL}chebi/terms/"))
        self.assertTrue(url.endswith("CHEBI_15377"))

    def test_http_error(self):
        url = f"{datasource.UNIPROT_BASE_URL}XXXXXX.json"
        with mock.patch("httpx2.get",
                        return_value=_response(url, 400, json={})):
            with self.assertRaises(httpx2.HTTPStatusError):
                fetch_uniprot_json("XXXXXX")

    def test_pubchem_assay_url(self):
        with mock.patch("httpx2.get", side_effect=lambda url, **_:
                        _response(url, json={"label": "x"})) as get:
            datasource.pubchem_assay("862")
        self.assertEqual(
            get.call_args.args[0],
            f"{datasource.PUG_BASE_URL}assay/aid/862/concise/CSV")


class TestUniProtTerms(unittest.TestCase):
    def test_uniprot_go_terms(self):
        terms, names = uniprot_go_terms(uniprot_json())
        self.assertEqual(terms, {
            "Function": ["GO:0003677"],
            "Process": ["GO:0006915"],
            "Component": ["GO:0005634"]
        })
        self.assertEqual(names["GO:0006915"], "apoptotic process")
        self.assertEqual(names["GO:0005634"], "nucleus")

    def test_uniprot_go_terms_unexpected_format(self):
        uj = uniprot_json()
        uj["uniProtKBCrossReferences"][0]["properties"][0]["key"] = "Other"
        with self.assertRaises(ValueError):
            uniprot_go_terms(uj)

    def test_uniprot_chebi_terms(self):
        self.assertEqual(uniprot_chebi_terms(uniprot_json()),
                         ["CHEBI:15377", "CHEBI:29105"])

    def test_chebi_name(self):
        self.assertEqual(chebi_name({"label": "water"}), "water")


class TestLoadTable(unittest.TestCase):
    def test_load_table(self):
        with tempfile.TemporaryDirectory() as d:
            Path(d, "a.csv").write_text("CID,value\n1,0.5\n2,1.5\n")
            spec = dict(sourceType="CSV", sourcePath="a.csv",
                        sampleIdColumn="CID", valueColumn="value")
            s = datasource.load_table(spec, Path(d))
        self.assertEqual(s.to_dict(), {1: 0.5, 2: 1.5})


if __name__ == "__main__":
    unittest.main()
