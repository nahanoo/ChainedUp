import pandas as pd


def wells(plate: str, col: str):
    rows = "ABCDEFGH"
    """Return wells like P1A03..P1H03 (col must be two digits as a string)."""
    return [f"{plate}{r}{col}" for r in rows]


# calibration
CAL_P1_SUCCINATE = wells("P1", "03")
CAL_P1_GLUCOSE = wells("P1", "06")
CAL_P2_SUCCINATE = wells("P2", "03")
CAL_P2_GLUCOSE = wells("P2", "06")

# P1
AT_SUCCINATE = wells("P1", "01")
AT_GLUCOSE = wells("P1", "02")
OA_SUCCINATE = wells("P1", "04")
OA_GLUCOSE = wells("P1", "05")
ATOA_SUCCINATE = wells("P1", "07")
ATOA_GLUCOSE = wells("P1", "08")

# P2
AT_SUCCINATE_GLUCOSE = wells("P2", "01")
ATC2_SUCCINATE_GLUCOSE = wells("P2", "02")
OA_SUCCINATE_GLUCOSE = wells("P2", "04")
OAC2_SUCCINATE_GLUCOSE = wells("P2", "05")
ATOA_SUCCINATE_GLUCOSE = wells("P2", "07")
ATOAC2_SUCCINATE_GLUCOSE = wells("P2", "08")
SUCCINATE_CONCENTRATIONS = [
    3.75000,
    0.93750,
    0.23438,
    0.05859,
    0.01465,
    0.00366,
    0.00091,
    0.00023,
]
GLUCOSE_CONCENTRATIONS = [
    2.50000,
    0.62500,
    0.15625,
    0.03906,
    0.00977,
    0.00244,
    0.00061,
    0.00015,
]

AT_TIMEPOINTS = [0, 23.15, 29.4, 47.06, 54.05, 71.8, 77.13, 96.3]

OA_TIMEPOINTS = [0, 22.3, 28.3, 46, 52.8, 69.94, 77.7, 94]

AT_OA_TIMEPOINTS = [0, 24.3, 29.3, 46.68, 54.124, 70.4, 78, 95]
