# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 14:42:08 2023

@author: dl923 / leadbot
"""

import os
import pandas as pd
import numpy as np


cazyfile='cazy_data/ssd/biblio/cazy_data\cazy_data.txt'
cazydb=pd.read_csv(cazyfile, sep='\t', encoding='utf-8', names=['Family', 'Taxa', 'Strain','Accession','Source'], header=None)

n_sample=50

subset_df_unfiltered=pd.DataFrame()

###restrict to ncbi only
c=0
tot=0
ncbidf=cazydb[cazydb['Source']=='ncbi']
for x in cazydb['Family'].unique():
    if not '_' in x:
        if x in ['AA9','AA10','AA11','AA13','AA14','AA15','AA16','AA17']:
            tempdf=ncbidf[ncbidf['Family']==x]
            if len(tempdf)>=n_sample:
                tempdf=tempdf.sample(n_sample)
                tot+=n_sample
            else:
                tot+=len(tempdf)
            subset_df_unfiltered=pd.concat([subset_df_unfiltered, tempdf])
            c+=1
            if c%50==0:
                print(str(c))
print(c)
print(tot)

subset_df_unfiltered.to_csv("CAZydb_ALL_LPMOs_n50_unfiltered.csv", sep=',', encoding='utf-8')

#%% Restrict to single accession occurances only
### Restrict to NCBI only
c = 0
tot = 0
ncbidf = cazydb[cazydb['Source'] == 'ncbi']

# Identify accessions that occur in multiple families
accession_family_counts = ncbidf.groupby('Accession')['Family'].nunique()
single_family_accessions = accession_family_counts[accession_family_counts == 1].index

# Filter to only accessions with a single family annotation
ncbidf = ncbidf[ncbidf['Accession'].isin(single_family_accessions)]

# Target LPMO families
target_families = ['AA9','AA10','AA11','AA13','AA14','AA15','AA16','AA17']

subset_df_filtered = pd.DataFrame()
available_counts = {}
selected_counts = {}

for x in cazydb['Family'].unique():
    if not '_' in x and x in target_families:
        tempdf = ncbidf[ncbidf['Family'] == x]
        available_count = len(tempdf)
        available_counts[x] = available_count

        if available_count >= n_sample:
            selected = tempdf.sample(n_sample)
            selected_counts[x] = n_sample
            tot += n_sample
        else:
            selected = tempdf
            selected_counts[x] = available_count
            tot += available_count

        subset_df_filtered = pd.concat([subset_df_filtered, selected])
        c += 1
        if c % n_sample == 0:
            print(str(c))

print(f"\nTotal families processed: {c}")
print(f"Total sequences collected: {tot}")

print("\nFamily accession summary (after filtering):")
for fam in target_families:
    available = available_counts.get(fam, 0)
    selected = selected_counts.get(fam, 0)
    print(f"{fam}: {available} available, {selected} added to subset")

# Save result
subset_df_filtered.to_csv("CAZydb_ALL_LPMOs_n50_filtered.csv", sep=',', encoding='utf-8')

#%% Download from NCBI
#################################################################
############### Download from NCBI using Entrez #################
#################################################################
from Bio import Entrez, SeqIO
from tqdm import tqdm  

outfile='CAZydb_ALL_LPMOs_n50_unfiltered.fasta'
source_dataframe=subset_df_unfiltered

def fetch_protein_sequences(accession_numbers, output_file):
    Entrez.email = "daniel.leadbeater@york.ac.uk"

    # Fetch protein sequences from NCBI
    handle = Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta", retmode="text")

    records = []
    for record in tqdm(SeqIO.parse(handle, "fasta"), desc="Downloading"):
        records.append(record)

    handle.close()

    # Write protein sequences to the output file
    SeqIO.write(records, output_file, "fasta")

    return records

def identify_missing_sequences(input_accession_numbers, downloaded_records):
    downloaded_accessions = set(record.id for record in downloaded_records)
    missing_accessions = set(input_accession_numbers) - downloaded_accessions

    return missing_accessions

if __name__ == "__main__":
    # List of accession numbers to download
    accession_numbers = list(source_dataframe[source_dataframe['Source']=='ncbi']['Accession'])
    # Output file name
    output_file = outfile

    # Fetch and save protein sequences, and get the count
    downloaded_records = fetch_protein_sequences(accession_numbers, output_file)
    sequence_count = len(downloaded_records)

    # Identify missing sequences
    missing_sequences = identify_missing_sequences(accession_numbers, downloaded_records)

    print(f"{sequence_count} protein sequences downloaded and saved to {output_file}")

    if missing_sequences:
        print("Missing sequences:")
        for accession in missing_sequences:
            print(accession)
    else:
        print("All requested sequences were downloaded.")
    
#%%Split multifasta into families
# Load your dataframe
from collections import defaultdict

target_subset_df = subset_df_filtered
target_fasta =  'CAZydb_ALL_LPMOs_n50_filtered.fasta'
outfolder="Filtered_N50"

if not outfolder in os.listdir():
    os.mkdir(outfolder)
# Create mapping of accession to family
accession_to_family = dict(zip(target_subset_df['Accession'], target_subset_df['Family']))

# Load all sequences from the multifasta file
fasta_sequences = list(SeqIO.parse(target_fasta, "fasta"))  # adjust filename as needed

# Prepare dictionary to group sequences by family
family_sequences = defaultdict(list)

# Process each sequence
for record in fasta_sequences:
    acc = record.id.split('|')[0]  # adjust if accessions have different format

    if acc in accession_to_family:
        fam = accession_to_family[acc]
        # Update description
        record.description = f"Family:{fam}_{record.description}"
        # Group by family
        family_sequences[fam].append(record)

# Write one FASTA file per family
for fam, records in family_sequences.items():
    filename = f"{fam}_sequences.fasta"
    SeqIO.write(records, os.path.join(outfolder,filename), "fasta")
    print(f"Wrote {len(records)} sequences to {filename}")


#%% Split multifasta

from Bio import SeqIO
import os
output_folder='Unfiltered_n50_per_family'
input_multifasta='CAZydb_ALL_LPMOs_n50_unfiltered.fasta'

def split_multifasta(input_file, output_folder):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read the multifasta file
    records = SeqIO.parse(input_file, "fasta")

    # Split and save each sequence into individual fasta files
    for record in records:
        sequence_id = record.id.replace("|", "_")  # Replace characters that are not allowed in filenames
        output_file = os.path.join(output_folder, f"{sequence_id}.fasta")

        with open(output_file, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":

    # Output folder for individual fasta files
    if not output_folder in os.listdir():
        os.mkdir(output_folder)
    # Run the split_multifasta function
    split_multifasta(input_multifasta, output_folder)

    print(f"Multifasta file '{input_multifasta}' has been split into individual fasta files in '{output_folder}'.")

