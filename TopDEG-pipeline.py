import pandas as pd
import numpy as np
import plotly.express as px

# Load DESeq2 results (example format)
# Ensure your CSV file has columns: gene, log2FoldChange, pvalue, padj
data = pd.read_csv("DEG_results.csv")

# Filter data for meaningful analysis (optional)
# Remove rows with NaN adjusted p-values
data = data.dropna(subset=["padj"])

# Add a significance column based on thresholds
log2_fc_threshold = 1.0
padj_threshold = 0.05
data["significance"] = "Not Significant"
data.loc[(data["padj"] < padj_threshold) & (data["log2FoldChange"] > log2_fc_threshold), "significance"] = "Upregulated"
data.loc[(data["padj"] < padj_threshold) & (data["log2FoldChange"] < -log2_fc_threshold), "significance"] = "Downregulated"

# Convert p-values to -log10(p-value)
data["-log10(pvalue)"] = -data["pvalue"].apply(lambda p: -1 if p <= 0 else np.log10(p))

# Interactive volcano plot
fig = px.scatter(
    data,
    x="log2FoldChange",
    y="-log10(pvalue)",
    color="significance",
    hover_data=["gene"],  # Display gene names on hover
    color_discrete_map={
        "Upregulated": "red",
        "Downregulated": "blue",
        "Not Significant": "gray",
    },
    title="Interactive Volcano Plot",
    labels={"log2FoldChange": "Log2 Fold Change", "-log10(pvalue)": "-Log10(p-value)"},
)

# Add threshold lines (optional)
fig.add_shape(
    type="line",
    x0=-log2_fc_threshold,
    x1=-log2_fc_threshold,
    y0=0,
    y1=data["-log10(pvalue)"].max(),
    line=dict(color="black", dash="dot"),
)
fig.add_shape(
    type="line",
    x0=log2_fc_threshold,
    x1=log2_fc_threshold,
    y0=0,
    y1=data["-log10(pvalue)"].max(),
    line=dict(color="black", dash="dot"),
)
fig.add_shape(
    type="line",
    x0=data["log2FoldChange"].min(),
    x1=data["log2FoldChange"].max(),
    y0=-np.log10(padj_threshold),
    y1=-np.log10(padj_threshold),
    line=dict(color="black", dash="dot"),
)

# Show the plot
fig.show()

# Extract significant genes (Upregulated or Downregulated)
significant_genes = data[data["significance"].isin(["Upregulated", "Downregulated"])]

# Save the list of significant genes to CSV
significant_genes.to_csv("significant_genes.csv", index=False)

# Display a message
print(f"List of significant genes saved to 'significant_genes.csv'.")

# Load your CSV file (replace with your actual file path)
csv_file_path = "DEG_results.csv"  # Replace with your file path
data = pd.read_csv(csv_file_path)

# Ensure the CSV has columns for Gene ID and Log Fold Change (modify column names if needed)
# Let's assume the CSV has columns 'Gene_ID' and 'Log_Fold_Change'
# If the column names are different, update them accordingly
data = data[['gene', 'log2FoldChange']]  # Adjust column names if necessary

# Rank the genes based on Log Fold Change (descending order for highest values first)
data['Rank'] = data['log2FoldChange'].rank(ascending=False, method='min')

# Sort the data by rank
ranked_data = data.sort_values(by='Rank', ascending=True)

# Print the top 100 genes
top_100_genes = ranked_data.head(100)

# Display the top 100 genes
print(top_100_genes)

# Optionally, you can save the top 100 genes to a new CSV file
top_100_genes.to_csv('top_100_genes.csv', index=False)


from Bio import Entrez

# Set your email address (required by NCBI)
Entrez.email = "your-email@example.com"

def search_pisum_sativum(locus_tag):
    try:
        # Search for the locus tag in the Gene database for Pisum sativum
        print(f"Searching for {locus_tag} in Pisum sativum...")
        handle = Entrez.esearch(db="gene", term=f"{locus_tag}[Gene] AND Pisum sativum[Organism]", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        if record["Count"] == "0":
            return f"No gene found for {locus_tag}"
        
        # Fetch the gene details using the gene ID
        gene_id = record["IdList"][0]
        print(f"Gene found. Gene ID: {gene_id}")
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        gene_record = Entrez.read(handle)
        handle.close()


 # Extract relevant gene information with checks for missing fields
        gene_info = {
            "Gene ID": gene_id,
            "Locus": locus_tag
        }
        
        # Extract Gene Description
        gene_info["Description"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_desc", "Uncharacterized")

        # Extract Protein Description
        gene_info["Protein Description"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_name", {}).get("Gene-ref_locus", "Uncharacterized protein")


# Functional Annotation (Gene Ontology)
        gene_info["Gene Ontology"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_function", "No functional annotation")


        return gene_info

    except Exception as e:
        return f"Error fetching data for {locus_tag}: {e}"

# Load the top 100 genes from the CSV
top_100_genes = pd.read_csv('top_100_genes.csv')

# Create an empty list to store results
gene_details = []


# Loop through the top 100 genes and fetch their details
for locus_tag in top_100_genes['gene']:
    gene_info = search_pisum_sativum(locus_tag)
    
    if isinstance(gene_info, dict):
        gene_details.append(gene_info)
    else:
        print(gene_info)  # Print any error message if no gene is found

# Convert the gene details list into a DataFrame
gene_details_df = pd.DataFrame(gene_details)

# Optionally, save the results to a new CSV file
gene_details_df.to_csv('top_100_gene_details.csv', index=False)

print("The gene details have been saved to 'top_100_gene_details.csv'.")





import requests
from bs4 import BeautifulSoup
import csv

def check_kegg_ko(gene_id):
    """
    Checks if a gene ID has a KO assigned in KEGG.
    """
    # Construct the KEGG URL for the gene
    url = f"https://www.genome.jp/dbget-bin/www_bget?{gene_id}"
    
    try:
        # Send GET request to KEGG website
        response = requests.get(url)
        
        # If the request was successful
        if response.status_code == 200:
            # Parse the HTML content of the page
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Look for KO in the page (this might need adjustment based on KEGG's HTML structure)
            ko_section = soup.find_all('b', string='KO:')
            
            if ko_section:
                # If KO is found, extract and return the KO value
                ko_value = ko_section[0].find_next('a').string.strip()
                return ko_value
            else:
                # If no KO is found, return "No KO assigned"
                return "No KO assigned"
        else:
            return f"Error: Unable to fetch data for {gene_id} (Status code: {response.status_code})"
    except Exception as e:
        return f"Error: {e}"

def check_genes_for_ko(gene_ids):
    """
    Checks a list of gene IDs to see if they have a KO assigned.
    """
    results = {}
    for gene_id in gene_ids:
        ko_info = check_kegg_ko(gene_id)
        results[gene_id] = ko_info
        print(f"Gene ID: {gene_id} - KO Info: {ko_info}")
    
    return results

# Step 1: Read gene IDs from top_100_genes.csv
df = pd.read_csv('top_100_genes.csv')  # Update the path if needed
# Assuming the gene IDs are in the 'gene' column (adjust the column name if different)
gene_ids = df['gene'].tolist()

# Step 2: Check KO for each gene ID from the CSV
results = check_genes_for_ko(gene_ids)

# Step 3: Optionally, save the results to a CSV file
with open('gene_ko_results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Gene ID", "KO Info"])
    for gene_id, ko_info in results.items():
        writer.writerow([gene_id, ko_info])

print("Results saved to gene_ko_results.csv")






