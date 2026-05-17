import pandas as pd
from Bio.PDB import PDBParser
import os

# Load distances data
csv_path = 'results/Final_Energies_and_Distances_True.csv'
df = pd.read_csv(csv_path)
df_clean = df.dropna(subset=['Closest_Cys_ID'])

# Aggregate hits per protein
counts = df_clean.groupby(['Protein', 'Closest_Cys_ID']).size().reset_index(name='Count')
counts = counts.sort_values(by=['Protein', 'Count'], ascending=[True, False])

docking_stats = {}
for prot in df_clean['Protein'].unique():
    prot_data = counts[counts['Protein'] == prot]
    stats_str = []
    for _, row in prot_data.iterrows():
        stats_str.append(f"Cys-{int(row['Closest_Cys_ID'])} ({row['Count']} hits)")
    docking_stats[prot] = " | ".join(stats_str)

# Extract full sequence Cys IDs
parser = PDBParser(QUIET=True)
all_sequence_cys = {}
proteins = sorted(df_clean['Protein'].unique())

for prot in proteins:
    pdb_path = os.path.join('data', prot)
    if not os.path.exists(pdb_path):
        pdb_path = os.path.join('data', 'proteins', prot)
        
    if os.path.exists(pdb_path):
        try:
            structure = parser.get_structure(prot, pdb_path)
            cys_list = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.resname == 'CYS':
                            cys_list.append(str(residue.id[1]))
                            
            all_sequence_cys[prot] = ", ".join(cys_list)
        except Exception as e:
            all_sequence_cys[prot] = "error"
    else:
        all_sequence_cys[prot] = "not found"

# Export summary table
final_rows = []
for prot in proteins:
    final_rows.append({
        'Protein': prot,
        'Sequence_All_Cys_IDs': all_sequence_cys.get(prot, ""),
        'Docking_Hits': docking_stats.get(prot, "")
    })

final_df = pd.DataFrame(final_rows)
output_file = 'results/Protein_Cys_Summary.csv'
final_df.to_csv(output_file, index=False)