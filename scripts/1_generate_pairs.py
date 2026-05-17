import os
import pandas as pd

# Define directories
prot_dir = 'data/proteins'
sub_dir = 'data/substrates'

# Select target proteins A1 to A7
target_proteins = [f"A{i}.pdb" for i in range(1, 8)]
proteins = [os.path.join(prot_dir, f) for f in target_proteins]

# Extract all sdf ligands
substrates = [os.path.join(sub_dir, f) for f in os.listdir(sub_dir) if f.endswith('.sdf')]

# Generate pair combinations
pairs = []
for p in proteins:
    for s in substrates:
        prot_name = os.path.basename(p).replace('.pdb', '')
        sub_name = os.path.basename(s).replace('.sdf', '')
        comp_name = f"{prot_name}_{sub_name}"
        
        pairs.append({
            'complex_name': comp_name,
            'protein_path': p,
            'ligand_description': s,
            'protein_sequence': ''
        })

df = pd.DataFrame(pairs)
df.to_csv('data/A1_A7_pairs.csv', index=False)