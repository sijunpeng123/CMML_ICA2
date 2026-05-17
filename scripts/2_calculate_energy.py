import os
import subprocess
import pandas as pd
import pymol

# Load pair lists to locate protein and ligand paths
df = pd.read_csv('data/A1_A7_pairs.csv')
results_dir = 'results/A1_A7_predictions'
final_results = []

# Initialize PyMOL API
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

for index, row in df.iterrows():
    comp_name = row['complex_name']
    prot_path = row['protein_path']
    ligand_name = os.path.basename(row['ligand_description'])
    prot_name = os.path.basename(prot_path)
    
    rank1_path = os.path.join(results_dir, comp_name, 'rank1.sdf')
    
    if not os.path.exists(rank1_path):
        continue
        
    temp_pdb = os.path.join(results_dir, comp_name, 'prodigy_temp.pdb')
    
    # Merge and rename chains using PyMOL
    pymol.cmd.reinitialize()
    pymol.cmd.load(prot_path, "prot")
    pymol.cmd.load(rank1_path, "lig")
    pymol.cmd.alter("prot", "chain='A'")
    pymol.cmd.alter("lig", "chain='L'")
    pymol.cmd.alter("lig", "resn='UNK'")
    pymol.cmd.create("complex", "prot or lig")
    pymol.cmd.save(temp_pdb, "complex")
    
    # Calculate binding free energy via PRODIGY-LIG
    cmd = ["prodigy_lig", "-i", temp_pdb, "-c", "A", "L:UNK"]
    energy = None
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in res.stdout.split('\n'):
            if line.strip().startswith('prodigy_temp'):
                parts = line.split()
                if len(parts) >= 2:
                    energy = float(parts[-1])
    except Exception as e:
        pass
        
    final_results.append({
        'Complex_Name': comp_name,
        'Protein': prot_name,
        'Ligand': ligand_name,
        'Delta_G (kcal/mol)': energy
    })

# Save the final energies report
df_results = pd.DataFrame(final_results)
df_results.to_csv('results/Final_Energies_Report.csv', index=False)