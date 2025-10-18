import os
from glob import glob
from datetime import datetime, timedelta
import pandas as pd


# --- Time helpers ---
def _parse_hhmmss(s: str) -> datetime:
    s = s.strip()
    fmt = "%H:%M:%S" if ":" in s else "%H%M%S"
    return datetime.strptime(s, fmt)


def excel_time_diff(t1: str, t0: str, wrap_midnight: bool = True) -> float:
    """Return (t1 - t0) as an Excel time fraction of a day."""
    dt1 = _parse_hhmmss(t1)
    dt0 = _parse_hhmmss(t0)
    delta = dt1 - dt0
    if wrap_midnight and delta.total_seconds() < 0:
        delta += timedelta(days=1)
    return delta.total_seconds() / 86400.0  # Excel fraction of a day


# --- File/plate parsing helpers ---
def get_run_index(fname: str) -> int:
    # Example name: "P3T194-20-09-25 093253_28.1C.csv"
    left = os.path.basename(fname).split("-")[0]  # "P3T194"
    return int(left[left.find("T") + 1 :])


def get_clock_hhmmss(fname: str) -> str:
    # Part after the space up to the underscore: "093253"
    return os.path.basename(fname).split(" ")[1].split("_")[0]


def list_wells():
    return [row + str(col) for row in "ABCDEFGH" for col in range(1, 13)]


def parse_od600_block(file_path: str) -> dict:
    """
    Parse the first block whose last token is 'OD600:600'.
    Return dict { 'A1': float, ..., 'H12': float }.
    """
    wells = list_wells()
    data = {}

    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    # Find header lines that start with ',' and include 'Wave Length'
    # The OD block has 9 lines: header + 8 row lines (A..H)
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(",") and "Wave Length" in line:
            # Candidate block: take next 8 lines for A..H
            candidate = lines[i : i + 9]
            if len(candidate) < 9:
                i += 1
                continue
            # Check last token of row A for OD600:600
            row_A = candidate[1]
            last_field = row_A.split(",")[-1].strip()
            if last_field.startswith("OD600:600"):
                # Parse rows A..H
                row_map = {r[0]: r for r in candidate[1:9]}  # 'A'..'H'
                for r in "ABCDEFGH":
                    parts = row_map[r].split(",")
                    # parts[0] is row letter, parts[1:13] are columns 1..12, parts[13] is label
                    nums = parts[1:13]
                    for c, val in enumerate(nums, start=1):
                        well = f"{r}{c}"
                        try:
                            data[well] = float(val)
                        except ValueError:
                            data[well] = float("nan")
                return data
            else:
                # Not the OD600 block, skip to next header
                i += 1
        else:
            i += 1

    raise ValueError(f"OD600:600 block not found in {file_path}")


# --- Main assembly ---
def build_plate_df(input_folder: str, plate_prefix: str) -> pd.DataFrame:
    files = sorted(
        glob(os.path.join(input_folder, f"{plate_prefix}*.csv")), key=get_run_index
    )
    if not files:
        raise FileNotFoundError(
            f"No files found for plate '{plate_prefix}' in {input_folder}"
        )

    # Base time is the clock time of the first file
    t0 = get_clock_hhmmss(files[0])

    # Prepare rows
    rows = []
    for f in files:
        t_curr = get_clock_hhmmss(f)
        t_frac = excel_time_diff(t_curr, t0, wrap_midnight=True)
        od_vals = parse_od600_block(f)
        row = {"Time": t_frac}
        row.update(od_vals)
        rows.append(row)

    # Order columns: Time, then wells in A1..H12 order
    wells = list_wells()
    df = pd.DataFrame(rows)
    df = df[["Time"] + wells].sort_values("Time", kind="stable").reset_index(drop=True)
    return df


# --------- Configure & run ----------
input_folder = "data/251018_succinate_glucose_plate_reader/raw"
plates = ["P1", "P2", "P3"]  # adjust as needed
output_xlsx = "plate_reader_OD600.xlsx"  # one file, a sheet per plate
output_folder = "data/251018_succinate_glucose_plate_reader/parsed"
os.makedirs(output_folder, exist_ok=True)
output_xlsx = os.path.join(output_folder, "plate_reader_OD600.xlsx")

with pd.ExcelWriter(output_xlsx, engine="xlsxwriter") as writer:
    for plate in plates:
        df_plate = build_plate_df(input_folder, plate)
        df_plate.to_excel(writer, sheet_name=plate, index=False)
