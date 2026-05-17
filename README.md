# CMML ICA2: Integrated Docking and Distance Analysis for Catalytic Cysteine Binding Sites

This repository contains the code, processed outputs, and figures for my CMML3 ICA2 mini-project.

The project investigates whether AI-predicted docking poses and binding-energy information are sufficient to identify biologically meaningful aldehyde-cysteine interactions, and whether three-dimensional distance analysis can improve the interpretation of predicted protein-ligand complexes.

## Project Overview

Cysteine residues are important nucleophilic sites in covalent drug design. In this project, aldehyde-containing ligands were docked against seven protein targets using DiffDock. The resulting protein-ligand complexes were analysed using PyMOL, PRODIGY-LIG, geometric distance calculation, and sequence-level cysteine validation.

The main aim was to determine whether docked aldehyde groups were positioned near catalytic cysteine residues, and whether different cysteine proximity patterns could be detected across protein targets.

## Repository Structure

```text
CMML-ICA2/
├── README.md
├── scripts/
│   ├── 1_generate_pairs.py
│   ├── 2_calculate_energy.py
│   ├── 3_calculate_distance.py
│   ├── 4_cys_statistics.py
│   └── 5_new_figures.py
├── figures/
│   ├── 1_workflow.png
│   ├── Cys_pattern.png
│   ├── Panel_B_Energy.png
│   ├── Panel_C_Patterns.png
│   └── Supp_Figure_S1.png
└── results/
    ├── Final_Energies_and_Distances_True.csv
    └── Protein_Cys_Summary.csv
