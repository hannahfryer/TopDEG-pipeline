## DEGExplorer: : Visualising and Annotating Top Genes

This repository contains scripts for analyzing differential gene expression (DEG) data, identifying significant genes, and annotating genes using external databases such as NCBI and KEGG. The scripts are designed to streamline the analysis and provide interactive visualization options, ranked gene lists, and functional annotations

## Repository Overview
Scripts:
1. interactive_volcano_plot.py: Processes DESeq2 results to identify significant genes based on adjusted p-values and log2 fold-change thresholds.
Saves significant genes to a CSV file and creates an interactive volcano plot for visualization.
DEG Analysis and Gene Information Retrieval

2. significant_gene_details.py: Processes DESeq2 results to rank genes by log2 fold change. Extracts the top 100 significant genes and queries the NCBI Gene database for Pisum sativum (pea) to retrieve gene details such as description, protein information, and functional annotations.

3. KEGG_retriever.py: Queries the KEGG database to check if gene IDs have KEGG Orthology (KO) assignments.
Reads gene IDs from a CSV file, retrieves KO information, and saves the results to a new CSV file.


Input Files:
DEG_results.csv: Contains DESeq2 results with columns such as gene, log2FoldChange, pvalue, and padj.
significant_genes.csv: A filtered list of significant genes generated by the first script.
top_100_gene_details.csv: A CSV file with the top 100 ranked genes, containing a column named Gene ID.

For a GitHub README, the **Requirements** section should be clear, structured, and concise. Here's a recommended structure:

---

## Requirements
- Python version 3.8 or higher

#### 1. Clone the Repository
```bash
git clone https://github.com/hannahfryer/DEGExplorer.git
cd DEGExplorer
```

#### 2. Install Dependencies
Use `requirements.txt` for a streamlined setup:
```bash
pip install -r requirements.txt
```

Alternatively, manually install the required libraries:
```bash
pip install pandas plotly numpy biopython beautifulsoup4 requests
```

---

This structure ensures users understand what they need and how to set it up efficiently. Each step is separated for clarity, and both automated and manual installation options are provided.

## Usage

### Input
If NOT using data supplied in this repository, ensure your input CSV file (DEG_results.csv) contains the following columns:
- gene (gene identifiers)
- lof2FoldChange (log2 fold changes)
- pvalue (p-values)
- padj (adjusted p-values)

### Running the Script
Run the script with the DEG results file as input:
```bash
python DEGExplorer.py --input DEG_results.csv
```

### Output
The following files will be generated:
- `volcano_plot.html`: Interactive volcano plot.
- `top_100_degs.csv`: Top 100 DEGs by log2 fold change.
- `top_100_gene_ids.csv`: Gene IDs and metadata from NCBI.
- `kegg_results.csv`: KEGG Orthology search results.

## Example

1. Run the pipeline:
   ```bash
   python TopDEG-pipeline.py 
   ```

2. Outputs:
   - `volcano_plot.html`
   - `top_100_degs.csv`
   - `top_100_gene_ids.csv`
   - `kegg_results.csv`



## Limitations
- The KEGG search may not return KO terms for all genes.
- Accuracy of gene ID retrieval depends on input data quality.
- Assumes properly formatted DEG input file.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
