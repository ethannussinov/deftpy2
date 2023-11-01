import pandas as pd
from crystal_analysis import Crystal, Defect

class KumagaiAnalysis:
    def __init__(self, data_path):
        self.data_path = data_path
        self.defects_data = self.load_defects_data()

    def load_defects_data(self):
        try:
            df = pd.read_csv(f"{self.data_path}/defects.csv")
            return df
        except Exception as e:
            print(f"Error loading defects data: {e}")
            return pd.DataFrame()

    def analyze_defects(self):
        if self.defects_data.empty:
            print("No defects data to analyze.")
            return

        results = []
        for index, row in self.defects_data.iterrows():
            defect = Defect(row['defect_type'], row['charge_state'], row['supercell_size'])
            crystal = Crystal()
            formation_energy = crystal.calculate_formation_energy(defect)
            results.append({
                'defect_type': defect.defect_type,
                'charge_state': defect.charge_state,
                'formation_energy': formation_energy
            })
            print(f"Formation Energy for {defect.defect_type} in charge state {defect.charge_state}: {formation_energy}")

        return pd.DataFrame(results)

if __name__ == "__main__":
    analyzer = KumagaiAnalysis("/path/to/your/data")
    result_df = analyzer.analyze_defects()
    if not result_df.empty:
        result_df.to_csv("/path/to/your/output/results.csv", index=False)
