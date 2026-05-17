import pandas as pd
from Bio.PDB import PDBParser
import numpy as np
import os

# Load docking energies report
df = pd.read_csv('results/Final_Energies_Report.csv')
results_dir = 'results/A1_A7_predictions'

min_distances = []
target_cys_ids = []
parser = PDBParser(QUIET=True)

# Process each complex to find minimum distance
for index, row in df.iterrows():
    comp_name = row['Complex_Name']
    pdb_path = os.path.join(results_dir, comp_name, 'prodigy_temp.pdb')
    
    if not os.path.exists(pdb_path):
        min_distances.append("")
        target_cys_ids.append("")
        continue

    try:
        structure = parser.get_structure(comp_name, pdb_path)
        model = structure[0]
        
        cys_sg_atoms = []
        ligand_c_atoms = []
        ligand_o_atoms = []
        
        for chain in model:
            for residue in chain:
                if residue.resname == 'CYS' and 'SG' in residue:
                    cys_sg_atoms.append((residue.id[1], residue['SG'].get_coord()))
                
                if residue.resname == 'UNK':
                    for atom in residue:
                        element = atom.element.strip().upper()
                        if not element:
                            element = ''.join([c for c in atom.name if c.isalpha()])[0].upper()
                        if element == 'C':
                            ligand_c_atoms.append(atom)
                        elif element == 'O':
                            ligand_o_atoms.append(atom)
                            
        if not cys_sg_atoms or not ligand_c_atoms or not ligand_o_atoms:
            min_distances.append("")
            target_cys_ids.append("")
            continue

        # Locate aldehyde carbon
        aldehyde_c_coord = None
        min_c_o_dist = float('inf')
        for c_atom in ligand_c_atoms:
            for o_atom in ligand_o_atoms:
                dist = np.linalg.norm(c_atom.get_coord() - o_atom.get_coord())
                if dist < min_c_o_dist:
                    min_c_o_dist = dist
                    if dist < 1.6: 
                        aldehyde_c_coord = c_atom.get_coord()
                        
        if aldehyde_c_coord is None:
            aldehyde_c_coord = min(ligand_c_atoms, key=lambda c: min([np.linalg.norm(c.get_coord() - o.get_coord()) for o in ligand_o_atoms])).get_coord()

        # Calculate distance to nearest Cys
        absolute_min_dist = float('inf')
        closest_cys_id = None
        for cys_id, sg_coord in cys_sg_atoms:
            dist = np.linalg.norm(sg_coord - aldehyde_c_coord)
            if dist < absolute_min_dist:
                absolute_min_dist = dist
                closest_cys_id = cys_id
        
        min_distances.append(round(absolute_min_dist, 3))
        target_cys_ids.append(closest_cys_id)
            
    except Exception as e:
        min_distances.append("")
        target_cys_ids.append("")

# Save results
df['Closest_Cys_ID'] = target_cys_ids
df['Distance_to_CHO (A)'] = min_distances
df.to_csv('results/Final_Energies_and_Distances_True.csv', index=False)