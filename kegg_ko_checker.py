# ---------------------------------------------------------------------------
# Script: KEGG KO Assignment Checker
# Description:
# - Queries KEGG database to check if gene IDs have KEGG Orthology (KO) assignments.
# - Processes gene IDs from a CSV file and retrieves KO information for each gene.
# - Saves the results to a new CSV file for further analysis.
# ---------------------------------------------------------------------------

# -------------------------- Library Imports --------------------------
# Ensure the necessary libraries are installed:

import requests  # For sending HTTP requests to the KEGG website
from bs4 import BeautifulSoup  # For parsing HTML content from KEGG pages
import pandas as pd  # For reading and processing data in structured format
import csv  # For writing results to a CSV file
import re  # For using regular expressions to identify KEGG Orthology (KO) patterns

# ------------------------ Define name of organism ------------------------
prefix = "psat"  # Replace "psat" with the desired prefix for the organism

# ------------------------ Function Definitions ------------------------
def check_kegg_ko(gene_id):
    """
    Query KEGG to check if a specific Gene ID has a KO (KEGG Orthology) assigned.

    Args:
        gene_id (str): The Gene ID to search for in KEGG.

    Returns:
        str: KO value(s) if found, "No KO assigned" if explicitly stated, or an error message.
    """
    # Construct the KEGG URL
    url = f"https://www.genome.jp/entry/{prefix}:{gene_id}"
    
    try:
        # Send a GET request to retrieve the KEGG page
        response = requests.get(url)
        
        if response.status_code == 200:
            # Parse the HTML content of the page
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Extract all the text content from the page
            page_text = soup.get_text(separator="\n", strip=True)

            # Use a regular expression to find KO entries (format: KXXXXX)
            ko_matches = re.findall(r'\bK\d{5}\b', page_text)

            if ko_matches:
                # Return all unique KO values as a comma-separated string
                return ", ".join(set(ko_matches)) 
            else:
                # Return a message if no KO is found
                return "KO not assigned"
        else:
            # Return an error message if the HTTP request fails
            return f"Error: Status code {response.status_code}"
    
    except Exception as e:
        return f"Error: {e}"

# --------------------- Main Script Execution ---------------------

#Load Gene IDs from CSV
# Assuming 'top_100_gene_details.csv' contains a column named "Gene ID"
input_csv = 'top_100_gene_details.csv'  # Path to your input file
df = pd.read_csv(input_csv)

# Extract the list of Gene IDs from the DataFrame
gene_ids = df['Gene ID'].tolist()

# Query KEGG for KO assignments
print("Checking KEGG KO assignments for Gene IDs...")

# Create an empty list to store results
results = []

# Loop through each Gene ID and query KEGG
for gene_id in gene_ids:
    # Use the check_kegg_ko function to get KO info for each Gene ID
    ko_info = check_kegg_ko(gene_id)

    # Append the results to the list as a dictionary
    results.append({"Gene ID": gene_id, "KO Info": ko_info})
    print(f"Gene ID: {gene_id} - KO Info: {ko_info}")  # Print progress to console

# Save Results to CSV
output_csv = 'gene_ko_results.csv'  # Path to save the results
results_df = pd.DataFrame(results)  # Convert results to DataFrame
results_df.to_csv(output_csv, index=False)  # Save to CSV

# Print a confirmation message upon successful completion
print(f"Results saved to {output_csv}")
