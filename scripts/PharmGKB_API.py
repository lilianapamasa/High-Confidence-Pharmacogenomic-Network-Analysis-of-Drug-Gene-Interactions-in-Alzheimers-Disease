import pandas as pd
import requests
import os
import sys
import json
import argparse
import time

output_dir = "/Users/lilianapamasa/Capstone"

#load significant DEGs
deg_df = pd.read_csv("/Users/lilianapamasa/Capstone/Filtered_LimmaResults.csv")

deg_list = deg_df["GeneSymbol"].unique().tolist()

#handle API requests with exponential backoff
def get_variant_drug_interactions(gene_symbol, max_retries=5, base_delay=2):
    """Fetch variant-drug interactions from PharmGKB API for a specific gene with rate limit handling"""
    url = f"https://api.pharmgkb.org/v1/data/variantAnnotation?location.genes.symbol={gene_symbol}"
    headers = {"Accept": "application/json"}
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers)
            
            if response.status_code == 429:
                wait_time = base_delay * (2 ** attempt)  #exponential backoff
                print(f"Rate limit exceeded. Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
                continue  #retry the request
            
            response.raise_for_status()  #raise for other errors (400, 500, etc.)
            
            data = response.json()
            if data.get("data") and isinstance(data["data"], list):
                return data["data"]
            return None
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {gene_symbol}: {e}")
            return None
    
    print(f"Max retries reached for {gene_symbol}. Skipping...")
    return None

results = [] #initialize a list to store results

#query each gene with controlled request rates
for gene in deg_list:
    interactions = get_variant_drug_interactions(gene)
    
    if interactions:
        for interaction in interactions:
            #extract variant
            variant = interaction.get("location", {}).get("rsid", "Unknown Variant")

            #extract all related drugs
            drugs = [chem.get("name", "Unknown Drug") for chem in interaction.get("relatedChemicals", [])]
            
            #store results properly
            for drug in drugs:
                results.append({
                    "Gene": gene,
                    "Variant": variant,
                    "Drug": drug
                })
    
    #add a short delay between requests to avoid 429 errors
    time.sleep(1)

results_df = pd.DataFrame(results) #convert results to DataFrame

results_df.head(10) #display results

output_dir = "/Users/lilianapamasa/Capstone"
output_file = f"{output_dir}/PharmGKB_Variant_Drug_Interactions.csv"

results_df.to_csv(output_file, index=False)
