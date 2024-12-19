# -----------------------------------------------
# Script: DEG Analysis and Interactive Volcano Plot
# Description: Processes DESeq2 results, identifies significant genes, 
# and visualizes results using an interactive volcano plot.
# -----------------------------------------------

# -------------------------- Library Imports --------------------------
# Install necessary libraries (if not already installed)
# pip install pandas plotly numpy 

import pandas as pd  # For data manipulation and analysis
import numpy as np  # For numerical computations and data processing
import plotly.express as px  # For creating interactive plots

# -------------------------- Data Loading ----------------------------

# Load DESeq2 results (example format)
# Ensure your CSV file has columns: gene, log2FoldChange, pvalue, padj
data = pd.read_csv("DEG_results.csv")

# -------------------------- Data Filtering --------------------------
# Filter data for meaningful analysis 
# Remove rows with NaN adjusted p-values
data = data.dropna(subset=["padj"])

# -------------------------- Significance Classification ------------
## Add a significance column based on thresholds
# Define thresholds for significance testing
log2_fc_threshold = 1.0 # Log2 fold-change threshold for up/downregulation
padj_threshold = 0.05 # Adjusted p-value threshold for statistical significance. Can change this if you want to make more stringent or less strict

# Creates a new column to classify gene expression significance
data["significance"] = "Not Significant"

# Label genes as "Upregulated" 
data.loc[(data["padj"] < padj_threshold) & (data["log2FoldChange"] > log2_fc_threshold), "significance"] = "Upregulated"
# Label genes as "Downregulated"  
data.loc[(data["padj"] < padj_threshold) & (data["log2FoldChange"] < -log2_fc_threshold), "significance"] = "Downregulated"

# -------------------------- Extract Significant Genes --------------
# Filter for upregulated and downregulated genes
significant_genes = data[data["significance"].isin(["Upregulated", "Downregulated"])]

# Save significant genes to a CSV file
output_significant_file = "significant_genes.csv"
significant_genes.to_csv(output_significant_file, index=False)
print(f"Significant genes saved to '{output_significant_file}'.")

# -------------------------- Prepare for Visualization --------------
# Convert p-values to -log10(p-value) for better visualisation on the y-axis
# Negative values are set to -1 to avoid issues with log(0) or negative p-values
data["-log10(pvalue)"] = -data["pvalue"].apply(lambda p: -1 if p <= 0 else np.log10(p))

# -------------------------- Volcano Plot Creation -------------------

# Create an interactive volvano plot to visualise gene expression results
fig = px.scatter(
    data,
    x="log2FoldChange", # x-axis: magnitude of log2 fold-change
    y="-log10(pvalue)", # y-axis: significance (-log10 of p-value)
    color="significance", # colour points by their significance category
    hover_data=["gene"],  # Display gene names when hover mouse over points
    color_discrete_map={  # map categories to specific colours
        "Upregulated": "red",
        "Downregulated": "blue",
        "Not Significant": "gray",
    },
    title="Interactive Volcano Plot", # title of plot
    labels={"log2FoldChange": "Log2 Fold Change", "-log10(pvalue)": "-Log10(p-value)"},
)

# -------------------------- Add Threshold Lines ---------------------
# Add vertical threshold lines for log2 fold-change significance
# These lines help identify upregulated and downregulated genes
fig.add_shape(
    type="line",
    x0=-log2_fc_threshold, # Negative fold-change threshold (downregulated)
    x1=-log2_fc_threshold,
    y0=0,
    y1=data["-log10(pvalue)"].max(),
    line=dict(color="black", dash="dot"),
)
fig.add_shape(
    type="line",
    x0=log2_fc_threshold, # Positive golf-change threshold (upregulated)
    x1=log2_fc_threshold,
    y0=0,
    y1=data["-log10(pvalue)"].max(),
    line=dict(color="black", dash="dot"),
)

# Add horizontal threshold line for p-value significance
# This line represents the adjusted p-value threshold (padj)
fig.add_shape(
    type="line",
    x0=data["log2FoldChange"].min(), # span the entire x-axis range
    x1=data["log2FoldChange"].max(), 
    y0=-np.log10(padj_threshold), # threshold for adjusted p-values
    y1=-np.log10(padj_threshold),
    line=dict(color="black", dash="dot"),
)

# Display the interactive plot
fig.show()
