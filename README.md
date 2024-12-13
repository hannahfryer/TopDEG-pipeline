# TopDEG-pipeline

# DEG Analysis and Functional Annotation Pipeline

This repository contains a Python-based pipeline designed to streamline the downstream analysis of differential gene expression (DEG) data. The script provides an interactive volcano plot, ranks DEGs based on log2 fold change, retrieves gene IDs and metadata from NCBI, and searches for KEGG Orthology (KO) terms for pathway analysis.

## Features

1. **Interactive Volcano Plot**
   - Visualize DEGs with an interactive volcano plot.
   - Adjust thresholds for significance directly in the interface.

2. **Top DEG Selection**
   - Automatically ranks DEGs based on log2 fold change.
   - Outputs the top 100 DEGs into a CSV file for further analysis.

3. **NCBI Gene ID Retrieval**
   - Queries NCBI to fetch gene IDs and metadata for the top DEGs.
   - Retrieves gene descriptions, genomic locations, and other relevant details.

4. **KEGG Orthology (KO) Search**
   - Searches KEGG database for KO terms using the retrieved gene IDs.
   - Outputs results for pathway analysis.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/hannahfryer/TopDEG-pipeline.git
   cd deg-analysis-pipeline
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Input
Prepare a CSV file containing DEG analysis results with the following columns:
- `gene_name`
- `log2_fold_change`
- `p_value`

### Running the Script
Run the script with the DEG results file as input:
```bash
python deg_analysis_pipeline.py --input <deg_results.csv>
```

### Output
The following files will be generated:
- `volcano_plot.html`: Interactive volcano plot.
- `top_100_degs.csv`: Top 100 DEGs by log2 fold change.
- `top_100_gene_ids.csv`: Gene IDs and metadata from NCBI.
- `kegg_results.csv`: KEGG Orthology search results.

## Workflow

1. **Prepare DEG Results**: Ensure input CSV contains `gene_name`, `log2_fold_change`, and `p_value` columns.
2. **Visualize with Volcano Plot**: Use the interactive HTML file to explore significant DEGs.
3. **Retrieve Gene IDs**: Generate a CSV with NCBI metadata for the top 100 DEGs.
4. **KEGG Search**: Use the KO terms from the output for pathway analysis.

## Example

1. Run the pipeline:
   ```bash
   python TopDEG-pipeline.py --input DEG_results.csv
   ```

2. Outputs:
   - `volcano_plot.html`
   - `top_100_degs.csv`
   - `top_100_gene_ids.csv`
   - `kegg_results.csv`

## Dependencies
- Python 3.8+
- pandas
- matplotlib
- plotly
- Biopython
- Requests

Install dependencies with:
```bash
pip install -r requirements.txt
```

## Limitations
- The KEGG search may not return KO terms for all genes.
- Accuracy of gene ID retrieval depends on input data quality.
- Assumes properly formatted DEG input file.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contributions
Contributions are welcome! Please fork the repository and submit a pull request.

## Contact
For questions or issues, please create an issue in the repository or email `your_email@example.com`.

