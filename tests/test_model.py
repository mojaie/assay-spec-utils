import unittest

from pydantic import ValidationError

from assay_spec_utils.model import (
    Assay, AssaySpec, DataSource, Readout, Target
)


class TestTarget(unittest.TestCase):
    def test_default_source_type(self):
        t = Target(accessionId="P04637").model_dump()
        self.assertEqual(t["sourceType"], "UniProt")
        self.assertIsNone(t["name"])

    def test_extra_field_forbidden(self):
        with self.assertRaises(ValidationError):
            Target(accessionId="P04637", unknown="x")


class TestReadout(unittest.TestCase):
    def test_defaults(self):
        r = Readout(readoutId="lum").model_dump()
        self.assertEqual(r["readoutMode"], "inhibition")
        self.assertEqual(r["readoutRange"], (0, 100))

    def test_invalid_readout_mode(self):
        with self.assertRaises(ValidationError):
            Readout(readoutId="lum", readoutMode="unknown")


class TestDataSource(unittest.TestCase):
    def test_screener_valid(self):
        d = DataSource(sessionId="AS-0001", layerIndices=[1, None])
        self.assertEqual(d.sourceType, "Screener")

    def test_screener_requires_session_and_layers(self):
        with self.assertRaises(ValidationError):
            DataSource(layerIndices=[1])
        with self.assertRaises(ValidationError):
            DataSource(sessionId="AS-0001")

    def test_csv_valid(self):
        d = DataSource(
            sourceType="CSV", sourcePath="a.csv", sampleIdColumn="CID",
            valueColumns=["value"])
        self.assertEqual(d.sourceType, "CSV")

    def test_csv_requires_fields(self):
        with self.assertRaises(ValidationError):
            DataSource(sourceType="CSV", sampleIdColumn="CID",
                       valueColumns=["value"])
        with self.assertRaises(ValidationError):
            DataSource(sourceType="CSV", sourcePath="a.csv",
                       valueColumns=["value"])
        with self.assertRaises(ValidationError):
            DataSource(sourceType="CSV", sourcePath="a.csv",
                       sampleIdColumn="CID")


class TestAssay(unittest.TestCase):
    def test_value_type_by_value(self):
        a = Assay(assayId="dose", valueType="Z-score").model_dump()
        self.assertEqual(a["valueType"], "Z-score")

    def test_invalid_value_type(self):
        with self.assertRaises(ValidationError):
            Assay(assayId="dose", valueType="IC50")


class TestAssaySpec(unittest.TestCase):
    def test_datasources_required(self):
        with self.assertRaises(ValidationError):
            AssaySpec(protocolId="p", assayId="a", readoutId="r")


if __name__ == "__main__":
    unittest.main()
