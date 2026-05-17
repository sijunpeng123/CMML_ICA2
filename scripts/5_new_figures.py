import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
sns.set_theme(style="ticks", context="talk")
plt.rcParams['font.family'] = 'sans-serif'

# Define file paths
csv_path = 'results/Final_Energies_and_Distances_True.csv'
output_dir = 'results'

# Load data
if not os.path.exists(csv_path):
    print(f"Error: {csv_path} not found.")
    exit()

df = pd.read_csv(csv_path)

# Clean data
df_clean = df.dropna(subset=['Closest_Cys_ID']).copy()
df_clean['Closest_Cys_ID'] = df_clean['Closest_Cys_ID'].astype(int)

# Define protein subgroups
subtype_vicinal = ['A1.pdb', 'A3.pdb', 'A4.pdb', 'A7.pdb']
subtype_single = ['A5.pdb', 'A6.pdb']

# Plot Panel B: Binding Free Energy Distribution
plt.figure(figsize=(7, 6))
sns.boxplot(x='Protein', y='Delta_G (kcal/mol)', data=df_clean, palette="Set3", width=0.6)
sns.stripplot(x='Protein', y='Delta_G (kcal/mol)', data=df_clean, color='black', alpha=0.3, size=4, jitter=True)
plt.axhline(y=-6.0, color='red', linestyle='--', linewidth=1.5, label='Threshold (-6.0)')
plt.ylabel('Delta G (kcal/mol)')
plt.xlabel('Protein Target')
plt.legend()
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Panel_B_Energy.png'), dpi=300)
plt.close()

# Plot Panel C: Cysteine Proximity Patterns
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

# Subtype I: Single-Cys Pattern
df_s1 = df_clean[df_clean['Protein'].isin(subtype_single)]
if not df_s1.empty:
    counts_s1 = df_s1.groupby(['Protein', 'Closest_Cys_ID']).size().reset_index(name='Hits')
    sns.barplot(x='Closest_Cys_ID', y='Hits', hue='Protein', data=counts_s1, ax=axes[0], palette="Blues_r")
    axes[0].set_title('Subtype I: Stable Single-Cys Pattern (A5, A6)')
    axes[0].set_ylabel('Hit Frequency')

# Subtype II: Vicinal Dithiol Alternating Pattern
df_s2 = df_clean[df_clean['Protein'].isin(subtype_vicinal)]
if not df_s2.empty:
    counts_s2 = df_s2.groupby(['Protein', 'Closest_Cys_ID']).size().reset_index(name='Hits')
    sns.barplot(x='Closest_Cys_ID', y='Hits', hue='Protein', data=counts_s2, ax=axes[1], palette="Reds_r")
    axes[1].set_title('Subtype II: Vicinal Dithiol Alternating Pattern (A1, A3, A4, A7)')
    axes[1].set_ylabel('Hit Frequency')
    axes[1].set_xlabel('Nearest Cysteine Residue ID')

sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Panel_C_Patterns.png'), dpi=300)
plt.close()

# Plot Supp Figure S1: Euclidean Distance Distribution
plt.figure(figsize=(8, 6))
sns.stripplot(x='Protein', y='Distance_to_CHO (A)', data=df_clean, jitter=True, alpha=0.6, size=5, palette="viridis")
plt.axhline(y=5.0, color='red', linestyle=':', linewidth=2, label='5.0 A Threshold')
plt.ylabel('Euclidean Distance (A)')
plt.xlabel('Protein Target')
plt.legend()
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Supp_Figure_S1.png'), dpi=300)
plt.close()

print("Execution complete. PNG files saved in results directory.")