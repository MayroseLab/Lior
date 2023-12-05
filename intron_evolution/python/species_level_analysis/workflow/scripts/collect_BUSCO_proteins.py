"""
Given a list of species, collect
orthologs into FASTA files
(one file per BUSCO)
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
from Bio import SeqIO

in_species_list = sys.argv[1]
in_dir = sys.argv[2]
out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# Read species list
with open(in_species_list) as f:
    species_list = [l.strip() for l in f.readlines()]

# Create DF with BUSCO names and SeqRecords

def binom_species(s):
    genus, species = s.split('_')[:2]
    return f'{genus}_{species}'

busco_prot = []
for sp in species_list:
    sp_binom = binom_species(sp)
    prot_fasta = os.path.join(in_dir, sp, 'prot.canon.fasta')
    prot_records = SeqIO.to_dict(SeqIO.parse(prot_fasta, 'fasta'))
    for rec_id in prot_records:
        prot_records[rec_id].id = f'{sp_binom}__{rec_id}'
        prot_records[rec_id].description = ''
    sp_busco_stats = os.path.join(in_dir, sp, 'BUSCO.stats')
    sp_busco_stats_df = pd.read_csv(sp_busco_stats, sep='\t', usecols=[1,5])
    sp_busco_stats_df['species'] = sp_binom
    sp_busco_stats_df['species_full'] = sp
    sp_busco_stats_df['prot_record'] = sp_busco_stats_df['canon_mRNA'].apply(lambda id: prot_records[id] if not pd.isnull(id) else np.nan)
    busco_prot.append(sp_busco_stats_df)
busco_prot_df = pd.concat(busco_prot)

# Extract proteins per BUSCO

MIN_SIZE = 5

def mad(s):
    s = s.dropna()
    return median_abs_deviation(s.values)

def mad_outliers(s):
    smed = s.median()
    smad = mad(s)
    lower = smed - 3*smad
    upper = smed + 3*smad
    return lower, upper

def create_BUSCO_fasta(busco_id):
    busco_id_df = busco_prot_df.query('Busco_id == @busco_id').dropna()
    busco_id_prot_records = busco_id_df['prot_record'].to_list()
    busco_id_prot_lens = pd.Series([len(rec) for rec in busco_id_prot_records])
    mad_min, mad_max = mad_outliers(busco_id_prot_lens)
    busco_id_prot_records_mad = [rec for rec in busco_id_prot_records if len(rec) > mad_min and len(rec) < mad_max]
    busco_id_fasta = os.path.join(out_dir, busco_id + '.faa')
    if len(busco_id_prot_records_mad) < MIN_SIZE:
        busco_id_prot_records_mad = []
    SeqIO.write(busco_id_prot_records_mad, busco_id_fasta, 'fasta')

for busco_id in busco_prot_df['Busco_id'].unique():
    create_BUSCO_fasta(busco_id)
