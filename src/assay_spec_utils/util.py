

def convert_conc_units(value: float, unit: str)-> tuple[float, str]:
    """Convert concentration values to their base units (M, g/L or ratio).
    """
    if value < 0:
        raise ValueError("concentration should be a positive value")
    mol_factor = dict(
        M=1, mM=1e-3, uM=1e-6, nM=1e-9, pM=1e-12
    )
    gl_factor = {
        "g/L": 1, "mg/L": 1e-3, "ug/L": 1e-6, "ng/L": 1e-9, "pg/L": 1e-12,
        "mg/mL": 1, "ug/mL": 1e-3, "ng/mL": 1e-6, "pg/mL": 1e-9,
        "ug/uL": 1, "ng/uL": 1e-3, "pg/uL": 1e-6
    }
    if unit in mol_factor.keys():
        return (value * mol_factor[unit], "M")
    elif unit in gl_factor():
        return (value * gl_factor[unit], "g/L")
    elif unit == r"v/v%":
        if value > 100 or value < 0:
            raise ValueError(r"v/v% should be in the range of 0-100")
        return (value * 0.01, "ratio")
    elif unit == "dilution":
        if value < 1:
            raise ValueError("dilution should be larger than 1")
        return (1.0 / value, "ratio")
    elif unit == "ratio":
        if value > 1:
            raise ValueError("ratio should be smaller than 1")
        return (value, unit)
    else:
        raise ValueError(f"invalid concentration unit notation: {unit}")


def convert_time_units(value: float, unit: str) -> tuple[int, str]:
    """Convert time values to their base units (seconds).
    """
    if value < 0:
        raise ValueError("time should be a positive value")
    time_factor = dict(
        sec=1, min=60, hour=3600, day=86400
    )
    if unit in time_factor.keys():
        return (value * time_factor[unit], "sec")
    else:
        raise ValueError(f"invalid time unit notation: {unit}")


def is_convertible_to_int(value) -> bool:
    """Check if the value is convertible to integer
    """
    try:
        int(value)
    except ValueError:
        return False
    return True


def is_convertible_to_float(value) -> bool:
    """Check if the value is convertible to float
    """
    try:
        float(value)
    except ValueError:
        return False
    return True