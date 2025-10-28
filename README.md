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

