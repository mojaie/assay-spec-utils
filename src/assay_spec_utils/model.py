
from enum import Enum
from typing import Optional, Dict, Tuple, Any

from pydantic import BaseModel, model_validator


class TargetSourceType(str, Enum):
    uniprot = 'UniProt'


class DataSourceType(str, Enum):
    screener = 'Screener'
    csv = 'CSV'


class ReadoutMode(str, Enum):
    inhibition = 'inhibition'
    activation = 'activation'


class ValueType(str, Enum):
    """Suggested value type. Data source APIs will be designed to consider this type and retrieve data in the appropriate format according to its application.
    """
    ac50 = 'AC50'
    zscore = 'Z-score'
    rzscore = 'RZ-score'
    percentage = 'percentage'
    category = 'category'  # nominal


class Target(BaseModel, extra='forbid', use_enum_values=True):
    sourceType: TargetSourceType = TargetSourceType.uniprot
    accessionId: str
    name: Optional[str] = None


class Attribute(BaseModel, extra='forbid'):
    attributeId: str
    name: Optional[str] = None
    description: Optional[str] = None
    active: Optional[bool] = True
    terms: list[Tuple[str, str]] = []  # termId, termName
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []  # category, name, value, unit


class Readout(BaseModel, extra='forbid', use_enum_values=True):
    readoutId: str
    readoutMode: ReadoutMode = ReadoutMode.inhibition
    readoutRange: Optional[Tuple[float,float]] = (0, 100)
    targets: list[Target] = []
    attributes: list[str] = []
    terms: list[Tuple[str, str]] = []
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []


class ProtocolTemplate(BaseModel, extra='forbid'):
    templateId: str
    name: Optional[str] = None
    description: Optional[str] = None
    active: Optional[bool] = True
    targets: list[Target] = []
    attributes: list[str] = []
    terms: list[Tuple[str, str]] = []
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []
    readouts: list[Readout] = []


class SubExperiment(BaseModel, extra='forbid'):
    subExpKey: str
    subExpValues: list[Optional[str]] = []


class DataSource(BaseModel, extra='forbid', use_enum_values=True):
    sourceType: DataSourceType = DataSourceType.screener
    sessionId: Optional[str] = None  # Screener
    layerIndices: list[Optional[int]] = []  # Screener
    subExperiments: list[SubExperiment] = []  # Screener, optional
    sourcePath: Optional[str] = None  # CSV, Excel
    sampleIdColumn: Optional[str] = None  # CSV, Excel
    valueColumns: list[Optional[str]] = []  # CSV, Excel
    sampleMapping: Optional[str] = None

    @model_validator(mode='after')
    def check_source_type(self):
        if (self.sourceType == "Screener"
                and (self.sessionId is None or not self.layerIndices)):
            raise ValueError(
                "sessionId and layerIndices should be specified")
        if (self.sourceType == "CSV" and (
                self.sourcePath is None or self.sampleIdColumn is None or not self.valueColumns)):
            raise ValueError(
                "sourcePath, sampleIdColumn and valueColumns should be specified")
        return self


class Assay(BaseModel, extra='forbid', use_enum_values=True):
    assayId: str
    valueType: ValueType = ValueType.ac50
    attributes: list[str] = []
    terms: list[Tuple[str, str]] = []
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []
    datasources: list[DataSource] = []


class AssayProject(BaseModel, extra='forbid'):
    assayProjectVersion: str = "1.0"
    meta: Dict[str, Any] = {}  # description


class AssayAttributes(BaseModel, extra='forbid'):
    assayAttributesVersion: str = "1.0"
    name: Optional[str] = None
    meta: Dict[str, Any] = {}  # description
    items: list[Attribute] = []  # TODO: attributeId duplication check


class AssayTemplates(BaseModel, extra='forbid'):
    assayTemplatesVersion: str = "1.0"
    name: Optional[str] = None
    meta: Dict[str, Any] = {}  # description
    items: list[ProtocolTemplate] = []  # TODO: templateId duplication check


class AssayProtocol(BaseModel, extra='forbid'):
    assayProtocolVersion: str = "1.0"
    meta: Dict[str, Any] = {}  # description
    templateId: Optional[str] = None
    targets: list[Target] = []
    attributes: list[str] = []
    terms: list[Tuple[str, str]] = []
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []
    readouts: list[Readout] = []  # TODO: readoutId duplication check
    assays: list[Assay] = []  # TODO: assayId duplication check


class AssaySpecSubExperiment(BaseModel, extra='forbid'):
    subExpKey: str
    subExpValue: str


class AssaySpecData(BaseModel, extra='forbid', use_enum_values=True):
    sourceType: DataSourceType = DataSourceType.screener
    sessionId: Optional[str] = None  # Screener
    layerIndex: Optional[int] = None  # Screener
    subExperiment: list[AssaySpecSubExperiment] = []  # Screener, optional
    sourcePath: Optional[str] = None  # CSV, Excel
    sampleIdColumn: Optional[str] = None  # CSV, Excel
    valueColumn: Optional[str] = None  # CSV, Excel
    sampleMapping: Optional[str] = None


class AssaySpec(BaseModel, extra='forbid', use_enum_values=True):
    """Assay specification for data analysis.

    Assays may have a datasetId like protocolId(_assayId)(_readoutId).
    Activity values will be retrieved according to Assay.datasources field.
    """
    protocolId: str  
    assayId: str
    readoutId: str
    readoutMode: ReadoutMode = ReadoutMode.inhibition
    readoutRange: Optional[Tuple[float,float]] = (0, 100)
    valueType: ValueType = ValueType.ac50
    targets: list[Target] = []
    terms: list[str] = []
    parameters: list[Tuple[str, str, Any, Optional[str]]] = []
    datasources: list[AssaySpecData]
