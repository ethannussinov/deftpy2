import pandas as pd
from glob import glob
from crystal_analysis import Crystal
from pymatgen.core import Structure

class WitmanDataAnalyzer:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data_frame = self.load_data()

    def load_data(self):
        csv_paths = sorted(glob(f"{self.data_path}/csvs/*.csv"))
        poscar_path = f"{self.data_path}/poscars"
        data_frames = []
        for csv_path in csv_paths:
            try:
                df = pd.read_csv(csv_path)
                filename = csv_path.split("/")[-1][:-4]
                df['filename'] = filename
                df['defectid'] = filename.split(".")[-1]
                df = df[df['defectname'] == "V_O"]
                df['poscar'] = df['filename'].apply(lambda x: open(f"{poscar_path}/{x}_POSCAR_wyck").read())
                df['structure'] = df['poscar'].apply(lambda x: Structure.from_str(x, fmt="poscar"))
                data_frames.append(df)
            except Exception as e:
                print(f"Error processing file {csv_path}: {e}")

        if not data_frames:
            return pd.DataFrame()
        return pd.concat(data_frames).reset_index(drop=True)

    def analyze_data(self):
        if self.data_frame.empty:
            print("No data to analyze.")
            return

        # Your data analysis logic goes here
        pass

if __name__ == "__main__":
    analyzer = WitmanDataAnalyzer("/path/to/your/data")
    analyzer.analyze_data()
