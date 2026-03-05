import pandas as pd


def parse_species(meta_path, data_path, filters):
    meta = pd.read_csv(meta_path)
    data = pd.read_csv(data_path)

    dfs = []
    for species_out, substrate_out, meta_species, carbon in filters:
        subset = meta[
            (meta["species"] == meta_species) & (meta["carbon_source"] == carbon)
        ]

        for i, lg in enumerate(subset["linegroup"]):
            time_col = f"{lg}_time"
            meas_col = f"{lg}_measurement"
            if time_col not in data.columns or meas_col not in data.columns:
                raise KeyError(f"Missing columns for linegroup {lg}")

            df = pd.DataFrame({"timepoint": data[time_col], "OD": data[meas_col]})
            df["species"] = species_out
            df["substrate"] = substrate_out
            df["replicate"] = f"{substrate_out.lower()}_{i}"
            df = df.dropna(subset=["timepoint", "OD"])
            dfs.append(df)

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def parse_at():
    return parse_species(
        "metaod/at/metadata.csv",
        "metaod/at/measurements.csv",
        filters=[
            ("At", "Glucose", "At pellet glucose", "Glucose"),
            ("At", "Succinate", "At pellet succinate", "Succinate"),
        ],
    )


def parse_oa():
    return parse_species(
        "metaod/oa/metadata.csv",
        "metaod/oa/measurements.csv",
        filters=[
            ("Oa", "Glucose", "Oa", "Glucose fresh medium"),
            ("Oa", "Succinate", "Oa", "Succinate fresh medium"),
        ],
    )
