import unittest
from deftpy.crystal_analysis import Crystal
from pathlib import Path

class TestTernaryVacancies(unittest.TestCase):
    def load_poscar_data(self, compound_name):
        # Define the path to the POSCAR file for the given ternary compound
        poscar_file_path = Path(__file__).parent / 'data' / 'test_files' / f'OQMD_{compound_name}_POSCAR.txt'
        
        # Read the POSCAR file and store its content
        with open(poscar_file_path, 'r') as poscar_file:
            poscar_data = poscar_file.read()
        
        return poscar_data

    def run_test(self, compound_name, element_symbol):
        # Load POSCAR data for the given ternary compound
        poscar_data = self.load_poscar_data(compound_name)

        # Test for Eb and Vr values for the specified element in the compound
        crystal = Crystal(poscar_string=poscar_data)
        eb_value = crystal.bond_dissociation_enthalpies[0][element_symbol]
        vr_value = crystal.reduction_potentials[0][element_symbol]

        # Add your assertions here to check if the values are as expected
        self.assertAlmostEqual(eb_value, expected_eb_value_for_element_in_compound)
        self.assertAlmostEqual(vr_value, expected_vr_value_for_element_in_compound)

    def test_eb_and_vr_for_CaAlO3_Ca(self):
        self.run_test("CaAlO3", "Ca")

    def test_eb_and_vr_for_CaAlO3_Al(self):
        self.run_test("CaAlO3", "Al")

    def test_eb_and_vr_for_CaAlO3_O(self):
        self.run_test("CaAlO3", "O")

    # Add similar test methods for other elements in different ternary compounds

if __name__ == '__main__':
    unittest.main()
