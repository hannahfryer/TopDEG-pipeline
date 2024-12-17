# ------------------- Ranking Genes Based on Log Fold Change -------------------

import pandas as pd

# Load CSV fule into a DataFrame
csv_file_path = "significant_genes.csv"  # Replace with your file path and/or adjust name of input file
data = pd.read_csv(csv_file_path)

# Select relevant columns for analysis
# Ensure the CSV contains columns for 'gene' and 'log2FoldChange' 
data = data[['gene', 'log2FoldChange']]  # modify column names here if your file uses different names

# Rank genes based on Log2 Fold Change
# Genes with higher log2FoldChange values will receive a higher rank (descending order)
data['Rank'] = data['log2FoldChange'].rank(ascending=False, method='min')

# Sort the data by rank to ensure proper order
ranked_data = data.sort_values(by='Rank', ascending=True)

# Extract the top 100 genes based on rank
top_100_genes = ranked_data.head(100) # Adjust number here if you want higher or lower amount of filtered genes

# Display the top 100 genes to the console
print("Top 100 genes ranked by Log2 Fold Change:")
print(top_100_genes)

# Save the top 100 genes to a new CSV file
# This file will contain the genes with the highest Log2 Fold change values
top_100_genes.to_csv('top_100_genes.csv', index=False)

# Print a confirmation message to the user
print("Top 100 genes saved to 'top_100_genes.csv")

# -------------------- NCBI Email Configuration --------------------
# Set your email address (required by NCBI to identify users)
Entrez.email = "your-email@example.com"

# -------------------- Function Definition --------------------
def search_pisum_sativum(locus_tag):
    """
    Search for gene details of a specific locus tag in the NCBI Gene database 
    for Pisum sativum (pea) and extract relevant gene information.
    
    Args:
        locus_tag (str): The locus tag of the gene to search for.
    
    Returns:
        dict: A dictionary containing Gene ID, description, protein details, and functional annotation.
        str: An error message if the gene is not found or an issue occurs.
    """
    try:
        # Step 1: Search for the locus tag in the Gene database
        print(f"Searching for {locus_tag} in Pisum sativum...")
        handle = Entrez.esearch(db="gene", term=f"{locus_tag}[Gene] AND Pisum sativum[Organism]", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        # Check if any results are found
        if record["Count"] == "0":
            return f"No gene found for {locus_tag}"
        
        # Step 2: Fetch the gene details using the first Gene ID
        gene_id = record["IdList"][0]
        print(f"Gene found. Gene ID: {gene_id}")
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        gene_record = Entrez.read(handle)
        handle.close()

        # Step 3: Extract relevant gene information with fallback values for missing fields
        gene_info = {
            "Gene ID": gene_id,        # The NCBI Gene ID
            "Locus": locus_tag         # The input locus tag
        }

        # Extract Gene Description (default: "Uncharacterized")
        gene_info["Description"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_desc", "Uncharacterized")
        
        # Extract Protein Description (default: "Uncharacterized protein")
        gene_info["Protein Description"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_name", {}).get("Gene-ref_locus", "Uncharacterized protein")

        # Extract Functional Annotation (Gene Ontology) if available
        gene_info["Gene Ontology"] = gene_record[0].get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_function", "No functional annotation")

        return gene_info  # Return the dictionary containing gene details
    
    except Exception as e:
        # Handle any errors that occur during the API request
        return f"Error fetching data for {locus_tag}: {e}"


# -------------------- Load Top 100 Genes --------------------
# Load the list of top 100 genes from a CSV file
# The file 'top_100_genes.csv' must contain a column named 'gene' with locus tags
top_100_genes = pd.read_csv('top_100_genes.csv')

# -------------------- Initialize Data Storage --------------------
# Create an empty list to store detailed gene information
gene_details = []

# -------------------- Fetch Gene Details --------------------
# Loop through each gene (locus tag) in the top 100 list
for gene_ID in top_100_genes['Gene ID']:
    # Fetch gene information using the 'search_pisum_sativum' function
    gene_info = search_pisum_sativum(locus_tag)
    
    # Check if valid gene information was returned
    if isinstance(gene_info, dict):
        gene_details.append(gene_info)  # Append valid results to the list
    else:
        print(gene_info)  # Print any error message for genes not found

# -------------------- Convert Results to DataFrame --------------------
# Convert the list of gene details (dictionaries) into a pandas DataFrame
gene_details_df = pd.DataFrame(gene_details)

# -------------------- Save Results to CSV --------------------
# Save the gene details DataFrame to a new CSV file for further analysis
gene_details_df.to_csv('top_100_gene_details.csv', index=False)

# Print confirmation message
print("The gene details have been saved to 'top_100_gene_details.csv'.")
