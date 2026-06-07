WSME / SASA Scripts

 Here I describe a small pipeline for converting PDB files into .npz, calculating residue SASA, running WSME/SASA fitting, and analyzing output.

Requirements:
numpy, scipy, pandas, matplotlib, openpyxl, biopython

pdbTOnpz.py also requires pyrosetta. It is recommended to run it on Linux/WSL/Server, because PyRosetta can be difficult to set up on Windows.


pdbTOnpz.py
Converts a PDB file into a compressed .npz file used by the model.
The output contains the energy tensor T, contact matrices C4, C6, C8, residue information, C-alpha coordinates, and gap information.
How to run the code in powershell/command etc.: 
python pdbTOnpz.py --pdb pdbs/1abc.pdb --out all_npz/1abc_full.npz --protein_only
Output:
.npz files

npzToSASA.py
Calculates residue-level SASA from a PDB file and checks whether the residue order matches the .npz file.
Outputs .csv and .npy SASA files.
Example:
python npzToSASA.py --npz all_npz/1abc_full.npz --pdb pdbs/1abc.pdb
If the PDB has multiple chains:
python npzToSASA.py --npz all_npz/1abc_full.npz --pdb pdbs/1abc.pdb --chain X
Output:
1abc_residue_sasa.csv
1abc_residue_sasa.npy

Ver1.py
Main fitting script.
It loads .npz files, SASA .npy files, and an Excel table with experimental folding ΔG values.
It optimizes model parameters and saves per-protein results plus restart parameters.
Example:
python Ver1.py --npz_dir all_npz --sasa_dir allSASA --excel "67 monomers revamped.xlsx" --cutoff C8 --min_seq_sep 3 --maxfev 150
Quick test with the smaller amount of proteins:
python Ver1.py --npz_dir all_npz --sasa_dir allSASA --excel "67 monomers revamped.xlsx" --cutoff C8 --min_seq_sep 3 --maxfev 20 --protein_limit 5
Continue the test from previous parameters:
python Ver1.py --npz_dir all_npz --sasa_dir allSASA --excel "67 monomers revamped.xlsx" --restart_params_file restart_params_C8_sep3_maxfev150.txt --maxfev 150
Output:
global_sasa_fit_perAA_results_C8_sep3_maxfev150.csv
restart_params_C8_sep3_maxfev150.txt

Ver1max and Ver2max just use a bit different optimization technique, all the above stays true.

analyze_wsme_results.py
Analyzes fitting results from CSV and generates plots and summary files.
Example:
python analyze_wsme_results.py --csv global_sasa_fit_perAA_results_C8_sep3_maxfev150.csv --npz_dir all_npz --cutoff C8 --min_seq_sep 3 --out_dir wsme_analysis
Output:
predicted_vs_experimental_dG.png, dG_vs_contacts.png, activation_barriers.png, per_protein_errors.png, summary.txt, results_with_contacts.csv

Typical Workflow
Convert PDB to NPZ
python pdbTOnpz.py --pdb pdbs/1abc.pdb --out all_npz/1abc_full.npz --protein_only

Calculate SASA
python npzToSASA.py --npz all_npz/1abc_full.npz --pdb pdbs/1abc.pdb

Run fitting
python Ver1.py --npz_dir all_npz --sasa_dir allSASA --excel "67 monomers revamped.xlsx" --cutoff C8 --min_seq_sep 3 --maxfev 150

Analyze results
python analyze_wsme_results.py --csv global_sasa_fit_perAA_results_C8_sep3_maxfev150.csv --npz_dir a
