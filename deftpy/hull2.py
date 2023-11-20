from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.entries.computed_entries import ComputedEntry
import pandas as pd

# Load the structure from the POSCAR file
poscar = Poscar.from_file("OQMD_CaTiO3_POSCAR.txt")
structure = poscar.structure

# Compute the total energy for the CaTiO3 structure (placeholder function)
# probably see if the value is attainable from the MP api first
def compute_total_energy(structure):
    # Implement the actual computation or retrieve the computed energy from MP
    return -20.0  # Placeholder energy value

# Create a ComputedEntry for CaTiO3
ca_tio3_energy = compute_total_energy(structure)
ca_tio3_entry = ComputedEntry(structure.composition, ca_tio3_energy)

# Load other entries from Eb and Vr datasets (placeholder function)
# You would replace this with actual data processing
def load_other_entries(eb_file, vr_file):
    # Load the Eb and Vr datasets
    eb_df = pd.read_csv(eb_file)
    vr_df = pd.read_csv(vr_file)

    # Convert Eb and Vr values to total energies (placeholder)
    # You need to implement the actual conversion logic
    entries = []
    for index, row in eb_df.iterrows():
        comp = Composition(row['compound'])
        energy = convert_to_total_energy(row['Eb'], row['Vr'])
        entry = PDEntry(comp, energy)
        entries.append(entry)
    return entries

# Assuming the convert_to_total_energy function exists and Eb.csv and Vr.csv are in the same directory
other_entries = load_other_entries("Eb.csv", "Vr.csv")

# Combine the CaTiO3 entry with the other entries
all_entries = other_entries + [ca_tio3_entry]

# Generate the phase diagram
phase_diagram = PhaseDiagram(all_entries)

# Calculate the Ehull for CaTiO3
ehull = phase_diagram.get_e_above_hull(ca_tio3_entry)

print(f"Ehull for CaTiO3 is {ehull} eV")
