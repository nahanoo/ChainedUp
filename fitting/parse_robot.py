#!/usr/bin/env python3
"""
Parse BioTek plate reader CSV files into tidy Excel workbooks (one sheet per plate),
handling multi-day runs by using full datetimes from filenames, and assuming each CSV
contains exactly THREE blocks in order:
    1) OD (ends with 'OD600:600')
    2) GFP (ends with '479,520' in the last two columns)
    3) mCherry (ends with '579,616' in the last two columns)

What you get in output_folder:
    - plate_reader_OD600.xlsx
    - plate_reader_GFP.xlsx
    - plate_reader_mCherry.xlsx

Each workbook has one sheet per plate (P1, P2, ...), columns:
    Time (Excel day fraction), A1..H12, Time_hours
"""

import os
import re
from glob import glob
from datetime import datetime
from typing import Dict, List, Tuple

import pandas as pd


# ----------------------------- Time & filename helpers -----------------------------

# Example filename: "P3T194-20-09-25 093253_28.1C.csv"
_DT_REGEX = re.compile(r"T\d+-(\d{2}-\d{2}-\d{2})\s+(\d{6})")


def get_file_dt(fname: str) -> datetime:
    """Extract datetime from filename like '...-20-09-25 093253_...csv' -> 2025-09-20 09:32:53."""
    base = os.path.basename(fname)
    m = _DT_REGEX.search(base)
    if not m:
        raise ValueError(f"Cannot parse date/time from filename: {base}")
    date_str, time_str = m.groups()
    return datetime.strptime(f"{date_str} {time_str}", "%d-%m-%y %H%M%S")


# ----------------------------- Plate parsing helpers ------------------------------


def list_wells() -> List[str]:
    """Return wells in A1..H12 order."""
    return [row + str(col) for row in "ABCDEFGH" for col in range(1, 13)]


def _iter_candidate_blocks(lines: List[str]):
    """
    Iterate over 96-well blocks:
      header line that starts with ',' and contains 'Wave Length' (or 'Wavelength'),
      followed by 8 lines (rows A..H).
    Yields the 8 row lines (strings) for each block.
    """
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(",") and ("Wave Length" in line or "Wavelength" in line):
            candidate = lines[i : i + 9]
            if len(candidate) >= 9:
                yield candidate[1:9]  # the 8 lines A..H
                i += 9
                continue
        i += 1


def _rows_to_well_dict(row_lines: List[str]) -> Dict[str, float]:
    """
    Convert 8 CSV row lines (A..H) into { 'A1': float, ..., 'H12': float }.
    Ignores any trailing label columns (e.g., 'OD600:600' or wavelength numbers).
    """
    wells = list_wells()
    data: Dict[str, float] = {}
    row_map = {r[0]: r for r in row_lines}  # keyed by 'A'..'H'
    for r in "ABCDEFGH":
        parts = row_map[r].split(",")
        # parts[0] is row letter; parts[1:13] are columns 1..12; any extras are labels
        nums = parts[1:13]
        for c, val in enumerate(nums, start=1):
            well = f"{r}{c}"
            try:
                data[well] = float(val)
            except ValueError:
                data[well] = float("nan")
    # Ensure all wells present
    for w in wells:
        data.setdefault(w, float("nan"))
    return data


def parse_three_blocks_in_order(
    file_path: str,
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, float]]:
    """
    Parse the FIRST THREE blocks found in the file, in order:
      1) OD, 2) GFP, 3) mCherry
    Returns three dicts (od_dict, gfp_dict, mcherry_dict).
    Raises if fewer than 3 blocks are present.
    """
    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    blocks = []
    for block_rows in _iter_candidate_blocks(lines):
        blocks.append(block_rows)
        if len(blocks) >= 3:
            break

    if len(blocks) < 3:
        raise ValueError(
            f"Expected at least 3 measurement blocks in {file_path}, found {len(blocks)}"
        )

    od = _rows_to_well_dict(blocks[0])
    gfp = _rows_to_well_dict(blocks[1])
    mch = _rows_to_well_dict(blocks[2])
    return od, gfp, mch


def build_plate_dfs(
    input_folder: str, plate_prefix: str
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build three DataFrames for one plate prefix (e.g., 'P3'), ordered by true datetime.
    Each DF has columns: 'Time' (Excel day fraction), A1..H12, 'Time_hours'.
    Returns: (df_od, df_gfp, df_mcherry).
    """
    pattern = os.path.join(input_folder, f"{plate_prefix}*.csv")
    files = sorted(glob(pattern), key=get_file_dt)
    if not files:
        raise FileNotFoundError(
            f"No files found for plate '{plate_prefix}' in '{input_folder}'. Looked for: {pattern}"
        )

    dt0 = get_file_dt(files[0])
    wells = list_wells()

    rows_od, rows_gfp, rows_mch = [], [], []

    for f in files:
        dt = get_file_dt(f)
        t_frac = (dt - dt0).total_seconds() / 86400.0  # Excel day fraction
        try:
            od_vals, gfp_vals, mch_vals = parse_three_blocks_in_order(f)
        except Exception as e:
            # Fail fast so you notice if a file is malformed
            raise RuntimeError(f"While parsing {f}: {e}") from e

        row_common = {"Time": t_frac}
        rows_od.append({**row_common, **od_vals})
        rows_gfp.append({**row_common, **gfp_vals})
        rows_mch.append({**row_common, **mch_vals})

    def _finish(rows: List[Dict[str, float]]) -> pd.DataFrame:
        df = pd.DataFrame(rows)
        df = (
            df[["Time"] + wells]
            .sort_values("Time", kind="stable")
            .reset_index(drop=True)
        )
        df["Time_hours"] = df["Time"] * 24.0
        return df

    return _finish(rows_od), _finish(rows_gfp), _finish(rows_mch)


# --------------------------------- Main script ------------------------------------

if __name__ == "__main__":
    # --------- Configure ----------
    input_folder = "data/251018_succinate_glucose_plate_reader/raw"
    plates = ["P1", "P2", "P3"]  # adjust as needed

    output_folder = "data/251018_succinate_glucose_plate_reader/parsed"
    os.makedirs(output_folder, exist_ok=True)

    # --------- Export per channel ----------
    od_xlsx = os.path.join(output_folder, "plate_reader_OD600.xlsx")
    gfp_xlsx = os.path.join(output_folder, "plate_reader_GFP.xlsx")
    mch_xlsx = os.path.join(output_folder, "plate_reader_mCherry.xlsx")

    with pd.ExcelWriter(od_xlsx, engine="xlsxwriter") as w_od, pd.ExcelWriter(
        gfp_xlsx, engine="xlsxwriter"
    ) as w_gfp, pd.ExcelWriter(mch_xlsx, engine="xlsxwriter") as w_mch:

        for plate in plates:
            df_od, df_gfp, df_mch = build_plate_dfs(input_folder, plate)
            df_od.to_excel(w_od, sheet_name=plate, index=False)
            df_gfp.to_excel(w_gfp, sheet_name=plate, index=False)
            df_mch.to_excel(w_mch, sheet_name=plate, index=False)

    print(f"✅ Wrote:\n  - {od_xlsx}\n  - {gfp_xlsx}\n  - {mch_xlsx}")
