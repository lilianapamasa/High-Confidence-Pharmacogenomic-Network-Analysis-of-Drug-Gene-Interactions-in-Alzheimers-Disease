#DGIdb GraphQL API
DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"

#initialize a list to store results
results = []

def get_drug_gene_interactions(gene):

    query = """
    query getInteractions($gene: [String!]!) {
      genes(names: $gene) {
        nodes {
          interactions {
            drug {
              name
              conceptId
            }
            interactionScore
            interactionTypes {
              type
              directionality
            }
            interactionAttributes {
              name
              value
            }
            publications {
              pmid
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }
    """
    variables = {"gene": [gene]}
    
    try:
        response = requests.post(
            DGIDB_GRAPHQL_URL,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"}
        )
        
        if response.status_code == 200:
            data = response.json()
            gene_nodes = data.get("data", {}).get("genes", {}).get("nodes", [])
            return gene_nodes[0].get("interactions", []) if gene_nodes else []
        else:
            print(f"Error {response.status_code}: {response.text}")
            return []
    
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return []

  #query each gene
for gene in deg_list:
    interactions = get_drug_gene_interactions(gene)
    
    if interactions:
        for interaction in interactions:
            #extract drug name and concept ID
            drug_name = interaction.get("drug", {}).get("name", "Unknown Drug")

            #extract interaction score
            score = interaction.get("interactionScore", "N/A")

            #extract interaction type
            interaction_types = [
                f"{itype['type']} ({itype.get('directionality', 'N/A')})"
                for itype in interaction.get("interactionTypes", [])
            ]
            interaction_type = ", ".join(interaction_types) if interaction_types else "Unknown"

            results.append({
                "Gene": gene,
                "Drug": drug_name,
                "Interaction Score": score,
                "Interaction Type": interaction_type
            })
    
    time.sleep(1)

#convert results to DataFrame
results_df = pd.DataFrame(results)

output_dir = "/Users/lilianapamasa/Capstone"
output_file = f"{output_dir}/DGIdb_gene_drug_interactions.csv"

results_df.to_csv(output_file, index=False)
