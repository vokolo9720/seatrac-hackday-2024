"""
python pre_pseudobulk.py \
  --sample_col biosample_id \
  --cluster_col CoarseClustering \
  --output_file coarse_pseudobulk_counts.csv

python pre_pseudobulk.py \
  --sample_col biosample_id \
  --cluster_col SubclusteringV2 \
  --output_file subclustering_pseudobulk_counts.csv


"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import sys
from fg_shared import _fg_data
from os.path import join as opj



def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Pseudo-bulk analysis from single-cell RNA-seq data.")
    parser.add_argument("--sample_col", default="sample_id", help="Column name in metadata for sample IDs")
    parser.add_argument("--cluster_col", default="cluster_id", help="Column name in metadata for cluster IDs")
    parser.add_argument("--output_file", default="pseudo_bulk_counts.csv", help="Output CSV file for pseudo-bulk counts")
    
    args = parser.parse_args()

    mcols = ['biosample_id', 'donor_id', 'array number', 'Sample Name',
               'Sample type', 'Time point of sampling',
               'Infusion before 2nd Mtb infection anti CD4 or IgG',
               'Naïve or Primary Infection or Reinfection at sample time', 'Group',
               'CFU Total ',
               'CoarseClustering',
               'SubclusteringV2', 'species']

    """Sum counts for each cell cluster within each sample for pseudo-bulk analysis."""
    base_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/bromley_etal')
    
    h5ad_fn = opj(base_fn, 'adata_Reinfection_CD4_Alexandria_SCP.h5ad')

    """Load single-cell data and metadata."""
    print("Loading single-cell data...")
    adata = sc.read_h5ad(h5ad_fn)

    # Group by sample and cell cluster
    grouped = adata.obs.groupby([args.sample_col, args.cluster_col])

    # Create a new DataFrame to store pseudo-bulk counts
    pseudo_bulk_counts = []
    gene_names = adata.var_names

    # Iterate over each group (sample, cluster) and sum counts
    for (sample_id, cluster_id), indices in grouped.indices.items():
        # Select the cells in this (sample, cluster) group
        sub_adata = adata[indices]
        
        # Sum the counts across cells
        counts_sum = np.asarray((np.exp(sub_adata.X.todense()) - 1).sum(axis=0)).flatten()
        
        # Store the sum in a dictionary
        pseudo_bulk_counts.append(pd.DataFrame({args.sample_col:[sample_id]*len(counts_sum),
                                    args.cluster_col:[cluster_id]*len(counts_sum),
                                    'counts':counts_sum,
                                    'gene':gene_names}))

    # Convert the pseudo-bulk dictionary to a DataFrame
    pseudo_bulk_df = pd.concat(pseudo_bulk_counts, axis=0)

    # Print summary and return
    print("Pseudo-bulk analysis completed!")
    print(pseudo_bulk_df.head())
    
    
    # pseudo_bulk_df = pd.merge(pseudo_bulk_df, adata.obs[mcols].drop_duplicates(), how='left', on=[args.sample_col, args.cluster_col])
    # Save the result
    pseudo_bulk_df.to_csv(args.output_file)
    adata.obs[mcols].drop_duplicates().to_csv(args.output_file.replace('.csv', '.meta.csv'))
    print(f"Pseudo-bulk counts saved to {args.output_file}")

if __name__ == "__main__":
    main()