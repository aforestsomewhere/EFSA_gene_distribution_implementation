import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
import numpy as np
import csv
import argparse
import sys


# Set up argument parser
parser = argparse.ArgumentParser(description="Process some input parameters.")
parser.add_argument('--output_dir', type=str, required=True, help='Output directory')
parser.add_argument('--genome_num', type=int, required=True, help='Number of genomes')
parser.add_argument('--species_name', type=str, required=True, help='Name of the species')
parser.add_argument('--blast_out', type=str, required=True, help='Output of blast')
parser.add_argument('--input_fasta', type=str, required=True, help='Input genes')
parser.add_argument('--identity_cutoff', type=int, required=True, help='Input genes')
parser.add_argument('--coverage_cutoff', type=int, required=True, help='Input genes')


# Parse arguments
args = parser.parse_args()
output_dir = args.output_dir
genome_num = args.genome_num
species_name = args.species_name
blast_out = args.blast_out
input_fasta = args.input_fasta
identity_cutoff = args.identity_cutoff
coverage_cutoff = args.coverage_cutoff


# Define path to BLAST output file
blast_output = os.path.join(os.getcwd(), blast_out)
amr_input = os.path.join(os.getcwd(), input_fasta)


# Read in the AMR_GENES_subt_fa.fasta file and create a dictionary of sequence lengths
seq_lengths = {}
for record in SeqIO.parse(amr_input, "fasta"):
    seq_lengths[record.id] = len(record.seq)

# Parse BLAST output file into a Pandas DataFrame
columns = ["query", "subject", "pct_identity", "alignment_length", "mismatches", "gap_opens",
           "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
blast_df = pd.read_csv(blast_output, sep="\t", header=None, names=columns)

# Add a new column to the blast_df DataFrame called query_length
blast_df["query_length"] = blast_df["query"].apply(lambda x: seq_lengths[x] if x in seq_lengths else len(x))

# Remove duplicates
blast_df = blast_df.drop_duplicates(subset=["query", "subject"])

#calculate the coverage
blast_df["coverage"] = (blast_df["alignment_length"] / blast_df["query_length"]) * 100

# Keeping only hits with coverage > 70 and Identity > 80
new_df = blast_df[(blast_df['coverage'] >= coverage_cutoff) & (blast_df['pct_identity'] >= identity_cutoff)]

if new_df.empty:
    print("No hits with enough coverage and identity, I will not generate any plot")
    sys.exit()

# Group the blast dataframe by query and calculate the percentage of strains that match each query
# the number of strains needs to be added manually
query_counts = blast_df.groupby("query")["subject"].nunique() / genome_num * 100
query_counts_thresh = new_df.groupby("query")["subject"].nunique() / genome_num * 100

# Save the DataFrames to  CSV file
blast_df.to_csv(os.path.join(os.getcwd(), output_dir, "blast_results.csv"), index=True)  # Set index=False to exclude the DataFrame index from the CSV file
new_df.to_csv(os.path.join(os.getcwd(), output_dir, "blast_results_filtered.csv"), index=True)  # Set index=False to exclude the DataFrame index from the CSV file


# Calculate the median for each query
result_df = query_counts.reset_index(name='pct_strains')
merged_df = pd.merge(result_df, blast_df, on="query", how="left")
median_identity_coverage = merged_df.groupby("query")[["pct_identity", "coverage"]].median().reset_index()
result_df = pd.merge(result_df, median_identity_coverage, on="query", how="left", suffixes=('', '_median'))
result_df.to_csv(os.path.join(os.getcwd(), output_dir, "catalogue_median.csv"), index=True)


result_df_tresh = query_counts_thresh.reset_index(name='pct_strains')
merged_df = pd.merge(result_df_tresh, blast_df, on="query", how="left")
median_identity_coverage = merged_df.groupby("query")[["pct_identity", "coverage"]].median().reset_index()
result_df_tresh = pd.merge(result_df_tresh, median_identity_coverage, on="query", how="left", suffixes=('', '_median'))
result_df_tresh.to_csv(os.path.join(os.getcwd(), output_dir, "catalogue_median_filtered.csv"), index=True)

# plot bar chart
fig, ax = plt.subplots(figsize=(20, 18))
query_counts.plot(kind="bar", ax=ax, color="#F47E13")
plt.xticks(rotation=90)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=2)
ax.set_xlabel("Gene")
ax.set_ylabel("Percentage of " + species_name + "strains with a match")
ax.set_title("Percentage of strains with a match for each gene")
plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), output_dir, "barchart_percentage_all.pdf"), format="pdf", dpi=400)

# plot bar chart
fig, ax = plt.subplots(figsize=(12, 11))
query_counts_thresh.plot(kind="bar", ax=ax, color="#44BDB7")
plt.xticks(rotation=90)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
ax.set_xlabel("Gene")
ax.set_ylabel("Percentage of " + species_name + "strains with a match")
ax.set_title("Percentage of strains with a match for each gene at 70% identity and 70% coverage thresholds")
plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), output_dir, "barchart_percentage_thresholds.pdf"), format="pdf", dpi=400)

 # Create heatmap using Seaborn % identity
heatmap_df = new_df.pivot(index="query", columns="subject", values="pct_identity")
plt.figure(figsize=(12, 11))
ax=sns.heatmap(heatmap_df, cmap="crest", annot=False, linewidth=0.5)
# set title
plt.title(species_name + ", Percentage identity, coverage higher than30%")
# set x and y axis labels
plt.xticks(rotation=90)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.set_xlabel("strains")
ax.set_ylabel("AMR gene")
plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), output_dir, "heatmap_identity_thresholds.pdf"), format="pdf", dpi=400)

 # Create heatmap using Seaborn % identity
heatmap_df = blast_df.pivot(index="query", columns="subject", values="pct_identity")
plt.figure(figsize=(20, 18))
ax=sns.heatmap(heatmap_df, cmap="crest", annot=False, linewidth=0.1)
# set title
plt.title(species_name + ", Percentage identity")
# set x and y axis labels
plt.xticks(rotation=90)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=6)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.set_xlabel("strains")
ax.set_ylabel("AMR gene")
plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), output_dir, "heatmap_identity_all.pdf"), format="pdf", dpi=400)

#scatter y= % of strains
g = sns.relplot(x='query', 
            y='pct_strains', 
            hue='pct_identity', 
            size='coverage',
            data=result_df,
            sizes=(40, 300), 
            alpha=.7, 
            palette='copper_r', 
            height=8, 
            aspect=8/8)

#added space to title; reduced x axis label size; rotated x axis labels.
title = "AMR genes in " + species_name
plt.title(title, fontsize=20, pad=20)
plt.subplots_adjust(top=0.85)
g._legend.texts[0].set_text('percentage of identity')
g._legend.set_bbox_to_anchor([1.01, .63])
plt.setp(g._legend.get_texts(), fontsize='6')
plt.xticks(rotation=90)
plt.savefig(os.path.join(os.getcwd(), output_dir, "scatter_all.pdf"), format="pdf", dpi=400)


plt.figure(figsize=(17, 10))
ax = sns.boxplot(x='query', y='pct_identity', hue="coverage", palette='winter' ,data=result_df)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.tick_params(axis='x', labelrotation=90, labelsize=6)
plt.savefig(os.path.join(os.getcwd(), output_dir, "box_DNA-all.pdf"), format="pdf", dpi=400)

plt.figure(figsize=(17, 10))
ax = sns.boxplot(x='query', y='pct_identity', hue="coverage", palette='winter' ,data=result_df_tresh)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.tick_params(axis='x', labelrotation=90, labelsize=6)
plt.savefig(os.path.join(os.getcwd(), output_dir, "box_DNA_thresholds.pdf"), format="pdf", dpi=400)
