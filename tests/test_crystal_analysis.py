import unittest
from deftpy.crystal_analysis import Crystal
from pathlib import Path

class TestOxygenVacancies(unittest.TestCase):
    def setUp(self):
        # Define the path to the POSCAR file
        test_file_path = Path(__file__).parent / 'data' / 'test_files' / 'OQMD_CaTiO3_POSCAR.txt'
        
        # Read the POSCAR file and store its content
        with open(test_file_path, 'r') as poscar_file:
            self.poscar_data = poscar_file.read()

    def test_eb_and_vr_for_CaTiO3_Ca(self):
        # Test for Eb and Vr values for Ca in CaTiO3
        crystal = Crystal(poscar_string=self.poscar_data)
        eb_value = crystal.bond_dissociation_enthalpies[0]['Ca']
        vr_value = crystal.reduction_potentials[0]['Ca']
        
        # Add your assertions here to check if the values are as expected
        self.assertAlmostEqual(eb_value, expected_eb_value_for_Ca_in_CaTiO3)
        self.assertAlmostEqual(vr_value, expected_vr_value_for_Ca_in_CaTiO3)

    def test_eb_and_vr_for_CaTiO3_Ti(self):
        # Test for Eb and Vr values for Ti in CaTiO3
        crystal = Crystal(poscar_string=self.poscar_data)
        eb_value = crystal.bond_dissociation_enthalpies[0]['Ti']
        vr_value = crystal.reduction_potentials[0]['Ti']
        
        # Add your assertions here to check if the values are as expected
        self.assertAlmostEqual(eb_value, expected_eb_value_for_Ti_in_CaTiO3)
        self.assertAlmostEqual(vr_value, expected_vr_value_for_Ti_in_CaTiO3)

    # Add similar test methods for other elements in your binary compound

if __name__ == '__main__':
    unittest.main()
