import pandas as pd


class Sample:
    def __init__(self, row):
        for col in row.index:
            setattr(self, col, row[col])


class Samples:
    def __init__(self):
        df = pd.read_excel("layout.ods", engine="odf", sheet_name="randomized_samples")
        self.samples = [Sample(df.iloc[i]) for i in range(len(df))]
