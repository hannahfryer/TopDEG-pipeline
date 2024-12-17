# ------------------------ Import Required Libraries ------------------------
import requests
from bs4 import BeautifulSoup
import pandas as pd
import csv

# ------------------------ Function Definitions ------------------------

import re  # Import regular expressions for pattern matching

def check_kegg_ko(gene_id):
    """
    Query KEGG to check if a specific Gene ID has a KO (KEGG Orthology) assigned.

    Args:
        gene_id (str): The Gene ID to search for in KEGG.

    Returns:
        str: KO value(s) if found, "No KO assigned" if explicitly stated, or an error message.
    """
    # Construct the KEGG URL
    url = f"https://www.genome.jp/entry/psat:{gene_id}"
    
    try:
        # Send GET request to KEGG
        response = requests.get(url)
        
        if response.status_code == 200:
            # Parse the HTML content
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Extract all the text content from the page
            page_text = soup.get_text(separator="\n", strip=True)

            # Check for explicit "No KO assigned" text
            if "no KO assigned" in page_text.lower():
                return "No KO assigned"

            # Use a regular expression to find KO entries (format: KXXXXX)
            ko_matches = re.findall(r'\bK\d{5}\b', page_text)

            if ko_matches:
                # Return all unique KO values as a comma-separated string
                return ", ".join(set(ko_matches))
            else:
                return "No KO assigned"
        else:
            return f"Error: Status code {response.status_code}"
    
    except Exception as e:
        return f"Error: {e}"

# ------------------------ Main Script Execution ------------------------

# Step 1: Load gene IDs from the input CSV file
# Assume the file 'top_100_genes.csv' has a column named 'gene' containing gene IDs
input_csv = 'top_100_genes.csv'  # Path to the input file
df = pd.read_csv(input_csv)
gene_ids = df['gene'].tolist()  # Extract the 'gene' column as a list of gene IDs

# Step 2: Check for KO assignments for each gene ID
print("Checking KEGG KO assignments for gene IDs...")
results = check_genes_for_ko(gene_ids)

# Step 3: Save the results to a new CSV file
output_csv = 'gene_ko_results.csv'  # Path to the output file
with open(output_csv, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Gene ID", "KO Info"])  # Write the header row
    
    # Write each gene ID and its corresponding KO information
    for gene_id, ko_info in results.items():
        writer.writerow([gene_id, ko_info])

# Print confirmation message
print(f"Results saved to {output_csv}")
