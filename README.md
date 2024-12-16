## DEGExplorer: : Visualising and Annotating Top Genes

This repository contains a Python-based pipeline designed for the processing, visualising, and annotating differential gene expression (DEG) data. The steps include generating volcano plots, ranking DEGs, retrieving
gene annotations, and identifying KEGG Orthology (KO) assignments. 

## Features

1. **Interactive Volcano Plot**:
   - Visualize DEGs with thresholds for log2 fold change and adjusted p-value.
   - Easily identify upregulated, downregulated, and non-significant genes.

2. **Gene Ranking**:
   - Rank genes by their log2 fold change and extract the top 100 genes for downstream analysis.

3. **Gene Annotation**:
   - Query NCBI databases to retrieve gene information and functional annotations for Pisum sativum (pea).

4. **KEGG KO Assignment Check**:
   - Identify KEGG Orthology assignments for genes using KEGG database queries.

5. **Output Generation**:
   - Save significant genes, ranked genes, and annotated gene information to CSV files for further analysis.
  
## Requirements
Install the following Python libraries before running the script: 
```bash
pip install pandas plotly numpy biopython beautifulsoup4 requests

## OR:

Clone the repository:
   ```bash
   git clone https://github.com/hannahfryer/DEGExplorer.git
   cd DEGExplorer
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

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

## Workflow

1. **Prepare DEG Results**: Ensure input CSV contains `gene_name`, `log2_fold_change`, and `p_value` columns.
2. **Visualize with Volcano Plot**: Use the interactive HTML file to explore significant DEGs.
3. **Retrieve Gene IDs**: Generate a CSV with NCBI metadata for the top 100 DEGs.
4. **KEGG Search**: Use the KO terms from the output for pathway analysis.

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
