
Assay-spec-utils
-------------------

This library is a Python utility for constructing structured, machine-readable assay specifications.

Installation
-----------------

- Requirements: Python 3.14+

```sh
pip install assay-spec-utils
```


Usage
-----------------

Note: The API is under development and may change.

### Specification files

Assay specifications are written as YAML files in the following directories (see `example/`).

```
src_dir/
├── attributes/  # Reusable sets of terms and parameters (e.g. reagents, microplates)
├── templates/   # Protocol templates (e.g. assay kits and detection methods)
└── protocols/   # Assay protocols. The file name is used as the protocolId.
```

Files whose names start with `_` are treated as drafts and ignored.

### Generating assay metadata

The `run` module processes the specification files in two steps.

```python
import logging
from pathlib import Path

from assay_spec_utils.run import fetch_target_info, process_all

logging.basicConfig(level=logging.INFO)

src_dir = Path("example")
dest_dir = Path("dest")

# Output directories must exist in advance.
for d in ["uniprot", "chebi", "processed"]:
    (dest_dir / d).mkdir(parents=True, exist_ok=True)

# 1. Download target information (UniProt and ChEBI JSON) used in the protocols.
#    Files already downloaded are skipped unless force_update=True.
fetch_target_info(src_dir, dest_dir)

# 2. Validate specifications and generate assay metadata (no network access).
process_all(src_dir, dest_dir)
```

`process_all` writes `dest_dir/processed/assays.json.gz`, which contains the following fields.

- `meta`: Creation timestamp.
- `assays`: Assay specifications for each combination of protocol, assay and readout, with resolved targets, terms, parameters and data sources.
- `targets`: GO terms (Function, Process and Component) and ChEBI terms for each target UniProt accession.
- `terms`: Dictionary of term IDs and names.

This library only records where the assay data are available (e.g. assay database sessions or CSV file paths) in `datasources`. Retrieving and validating the actual data are left to downstream applications.


License
-----------------

[MIT license](http://opensource.org/licenses/MIT)


Copyright
-----------------

(C) 2024-2026 Seiji Matsuoka