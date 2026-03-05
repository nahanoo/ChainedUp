import pandas as pd


def sort_raw():
    df = pd.read_csv("raw.csv", index_col="Name")
    meta = pd.read_excel(
        "2602246_sample_sheet_eric_ulrich_corrected.xlsx",
        sheet_name="randomized_samples",
    )
    order = meta["machine_filename"]
    out = pd.DataFrame(columns=df.columns)
    for i in order:
        row = df.loc[i]
        out.loc[i] = row
    out.to_csv("sorted_raw.csv")


def normalize_data():
    df = pd.read_excel("data.xlsx")
    df["glucose_area"] = df["glucose1_area"].fillna(0) + df["glucose2_area"].fillna(0)
    df["glucose_area_normalized"] = df["glucose_area"] / df["sorbitol_area"]
    df["succinate_area_normalized"] = df["succinate_area"] / df["norvaline_area"]
    df["ribose_area_normalized"] = df["ribose_area"] / df["sorbitol_area"]
    df.to_csv("normalized_data.csv", index=False)


normalize_data()
