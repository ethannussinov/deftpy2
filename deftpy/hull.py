from pymatgen.core import Structure
from pysipfenn import Calculator
import pandas as pd

eb_data = pd.read_csv('../data/Eb.csv')
vr_data = pd.read_csv('../data/Vr.csv')

def create_structure_from_poscar(poscar_path):
    with open(poscar_path, 'r') as file:
        poscar_content = file.read()
    structure = Structure.from_str(poscar_content, fmt="poscar")
    return structure

# Path to  POSCAR files
poscar_files = [
    '../data/OQMD_CaTiO3_POSCAR.txt',
    '../data/OQMD_CeO2_POSCAR.txt'
]

# Convert POSCAR files to pymatgen Structure objects
structures = [create_structure_from_poscar(path) for path in poscar_files]

# Initialize the pySIPFENN Calculator
calculator = Calculator()

descriptors = [calculator.calculate_KS2022([structure])[0] for structure in structures]

calculator.loadModels()

compatible_models = calculator.findCompatibleModels('KS2022')

predictions = [calculator.makePredictions(calculator.loadedModels, compatible_models, [descriptor])[0] for descriptor in descriptors]


ehull_values = [pred['Ehull'] for pred in predictions]

# Save the 'Ehull' predictions to a CSV file
ehull_df = pd.DataFrame({'Ehull': ehull_values})
ehull_df.to_csv('../data/Ehull_predictions.csv', index=False)

# Optionally, you can return the 'Ehull' values if you want to use them elsewhere in your code
def get_ehull_values():
    return ehull_values

if __name__ == "__main__":
    # Run the main process to generate 'Ehull' predictions
    get_ehull_values()
