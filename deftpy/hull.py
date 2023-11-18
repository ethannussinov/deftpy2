from pysipfenn import Calculator
import pysipfenn.descriptorDefinitions as desc_defs
from pymatgen.core import Structure
import pandas as pd

def calculate_ehull_for_structure(file_path):
    # Initialize the pySIPFENN Calculator
    calculator = Calculator()

    # Load the structure from the file
    with open(file_path, 'r') as file:
        structure = Structure.from_str(file.read(), fmt="poscar")

    # Generate the KS2022 descriptor for the structure
    # The descriptor generation function is called directly from the descriptorDefinitions module
    descriptor = desc_defs.KS2022.generate_descriptor(structure)

    # Load models compatible with 'KS2022' descriptors
    calculator.loadModels()
    print("Loaded models:", calculator.loadedModels.keys())


    # Find compatible models for the descriptor
    compatible_models = calculator.findCompatibleModels('KS2022')

    # Before making predictions
    print("Compatible models found:", compatible_models)
    print("Making predictions using loaded models and calculated descriptor...")

    # Make sure you're only using models that have been loaded successfully
    loaded_and_compatible_models = {model_name: calculator.loadedModels[model_name] 
                                    for model_name in compatible_models 
                                    if model_name in calculator.loadedModels}

    # Now use the filtered dictionary to make predictions
    predictions = calculator.makePredictions(loaded_and_compatible_models, list(loaded_and_compatible_models.keys()), [descriptor])

    # Initialize an empty list for ehull values
    ehull_values = []

    # Try to extract 'Ehull' values from predictions
    try:
        # Assuming that each prediction is a dictionary with 'Ehull' as one of the keys
        ehull_values = [pred['Ehull'] for pred in predictions]
    except KeyError:
        # If 'Ehull' key is not present in the predictions dictionary
        print("The 'Ehull' key was not found in the prediction output.")
    except TypeError:
        # If the predictions are not in an expected format (not a list or not a dictionary)
        print("The predictions were not in an expected format.")

    # Create a DataFrame with the 'Ehull' predictions
    ehull_df = pd.DataFrame({'Ehull': ehull_values})

    return ehull_df

# test the function with an example POSCAR file.
if __name__ == "__main__":
    example_file_path = '../data/test_files/OQMD_CaTiO3_POSCAR.txt'  # Replace with the actual file path
    ehull_df = calculate_ehull_for_structure(example_file_path)
    print(ehull_df)
