\# Protein Resetability (R) Analysis Pipeline



\[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



The `run\_protein\_R.py` script implements the \*\*Resetability (R)\*\* metric, a fast, geometry-based tool for analyzing protein backbone flexibility directly from single PDB/mmCIF structures.



This metric is designed to quickly identify dynamic regions, hinge points, and changes in structural rigidity associated with biological function, allostery, and drug binding.



\## Key Features



\* \*\*Metric Validation:\*\* R is benchmarked on classical systems like \*\*Hemoglobin T↔R\*\* allostery and \*\*PD-1\*\* antibody-induced rigidification.

\* \*\*Immunology Focus:\*\* Provides comparative analysis tools useful for studying \*\*HLA\*\* allele mechanics and \*\*TCR-HLA\*\* complex formation.

\* \*\*Fully Automated Output:\*\* Generates CSV summaries, R profile plots, Delta R ($\\Delta R$) comparison plots, and a draft manuscript (`.docx`).

\* \*\*No Local Setup Required:\*\* The entire pipeline can be run instantly via Google Colab.



\## Requirements



\* Python 3.7+

\* The necessary libraries are listed in `requirements.txt` (including Biopython, NumPy, Pandas, Matplotlib, and python-docx).



\## Setup \& Usage (Local Machine)

🧩 Required Benchmark Pairs

To fully reproduce the validation tests described in the manuscript and figures, the following protein structure pairs are required.
Each pair represents a biologically relevant conformational change, used to test the Resetability (R) metric.

Label	PDB IDs	Description	Biological Purpose
Hemoglobin_T_vs_R	2HBB (T state), 1AJ9 (R state)	Human hemoglobin in tense (T) and relaxed (R) forms	Classic allosteric transition test for ΔR sensitivity
PD1_Apo_vs_Bound	3RRQ, 3B71 (apo) → 5GGS (bound)	PD-1 receptor, free vs antibody-bound	Measures ligand-induced rigidification
HLA_A0201_vs_B27	1A1N (HLA-A0201), 1HSA (HLA-B27)	Two HLA alleles	Compares immune allotype structural dynamics
HLA_A2_TCR_vs_free	1A1N (free), 1AO7 (TCR-bound)	HLA–TCR complex formation	Captures receptor-induced rotation and stability shift

These structures are automatically loaded when the script detects their .cif or .mmCIF files under:

Protein_R/inputs/


For best results, include the following files in that folder:

1AJ9.cif
2HBB.cif
3RRQ.cif
3B71.cif
5GGS.cif
1A1N.cif
1HSA.cif
1AO7.cif


All of these can be downloaded directly from the RCSB Protein Data Bank
 by searching for their PDB IDs.
The script will automatically parse whichever subset is available and skip missing entries gracefully.

1\.  \*\*Install Dependencies:\*\*

&nbsp;   ```bash

&nbsp;   pip install -r requirements.txt

&nbsp;   ```

2\.  \*\*Prepare Input:\*\* Create a folder named `Protein\_R/inputs/` and place your protein structure files (`.pdb`, `.cif`, or `.cif.gz`) inside it.

3\.  \*\*Run Analysis:\*\*

&nbsp;   ```bash

&nbsp;   python run\_protein\_R.py

&nbsp;   ```

4\.  \*\*Check Output:\*\* All results (CSV files, PNG plots, DOCX report) are saved to the `Protein\_R/results/` directory.



\## 🚀 Instant Usage (Google Colab)



To run the full pipeline without installing anything locally, use the dedicated Colab notebook:



🔗 \*\*Google Colab Notebook:\*\* `\[https://colab.research.google.com/drive/1PZgQSTojclud9jZ6t_uKRNAFPeboC1lY?usp=sharing]`



The notebook will handle file uploads, package installation, and automatically download the final results zip file containing the plots and data.



\## Customization



The core analysis can be modified by editing the constants near the top of `run\_protein\_R.py`:



\* \*\*`WINDOW\_SIZE` / `WINDOW\_STEP`\*\*: Change the length and step-size of the sliding window analysis.

\* \*\*`PAIRS` Dictionary\*\*: \*\*Crucially\*\*, edit this dictionary to define which PDB IDs (or filename substrings) you want the script to compare (e.g., comparing your 'apo' structure ID with your 'bound' structure ID).



\## Citation



\*(Placeholder: Add how you would like others to cite your work—likely linking to the bioRxiv preprint and/or the GitHub repository.)\*



\## License



This project is licensed under the MIT License. See the \[LICENSE](LICENSE) file for details.

