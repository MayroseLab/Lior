"""
This is the final step of the
discepancies analysis pipeline,
that parses all pipeline outputs
and assigns one of several
possible causes to each discrepancy
"""

import argparse
import pandas as pd
import os
from pysam import VariantFile
from itertools import chain

def parse_psl(psl):
  """
  Parse a BLAT output file
  in psl format.
  """
  colnames = ["matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts"]
  df = pd.read_csv(psl, sep='\t', names=colnames, skiprows=5)
  return df

def create_data_dict(in_dir, ext, func):
  """
  Given a dir containing files
  <name>.<ext>, parses each file using
  func and returns a dict:
  {name1 : func_res1, name2 : func_res2, ...}
  """
  data_dict = {}
  for d in [f for f in os.listdir(in_dir) if f.endswith('.'+ext)]:
    sample = d.split('.')[0]
    data_dict[sample] = func(os.path.join(in_dir, d))
  return data_dict

def query_in_psl(psl_df, query, min_frac):
  """
  Returns a Series with the best
  hit of a given query.
  Query name must  be in qName column,
  and matches/qSize must be larger
  than min_frac. If no record passes
  the threshold - return empty, if
  multiple hits - return the one
  with highest matches/qSize
  """
  hits = psl_df.loc[(psl_df['qName'] == query) & (psl_df['matches']/psl_df['qSize'] > min_frac)]
  if hits.empty:
    return hits
  return hits.iloc[(hits['matches']/hits['qSize']).argmax()]

def parse_snpEff_bcf(bcf):
  """
  Parses a BCF file output of
  SnpEff to create a query-able
  data structure:
  {gene1: [{variant1_ANN}, {variant2_ANN}], gene2: [...]}
  """
  res = {}
  ann_keys = ("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos_cDNA.length","CDS.pos_CDS.length","AA.pos_AA.length","Distance","ERRORS_WARNINGS_INFO")
  for rec in VariantFile(bcf):
    for ann in rec.info['ANN']:
      ann_dict = dict(zip(ann_keys, ann.split('|')))
      gene = ann_dict["Feature_ID"]
      if gene not in res:
        res[gene] = []
      res[gene].append(ann_dict)
  return res

def is_pseudogene(gene, variant_annotation_dict):
  """
  Search for HIGH impact variants
  in a given gene, within a variant
  annotation dict.
  Returns True/False
  """
  if gene not in variant_annotation_dict:
    return False
  for ann in variant_annotation_dict[gene]:
    if ann['Annotation_Impact'] == 'HIGH':
      return True
  return False

def parse_og_tsv(og_tsv):
  """
  Read an Orthogroups TSV and create
  a dictionary:
  {gene1: OG, gene2: OG,...}
  """
  og_dict = {}
  with open(og_tsv) as f:
    f.readline()	# skip header
    for line in f:
      line = line.strip()
      fields = line.split('\t')
      og = fields[0]
      genes = list(chain.from_iterable([s.split(', ') for s in fields[1:]]))
      for gene in genes:
        og_dict[gene] = og
  return og_dict

def are_clustered_differently(gene1, gene2, og_orig_dict, og_mwop_dict):
  """
  Are gene1 and gene2 in the same OG
  in og_orig_dict and different OGs
  in og_mwop_dict? Return True/False
  """
  if og_orig_dict[gene1] == og_orig_dict[gene2] and og_mwop_dict[gene1] != og_mwop_dict[gene2]:
    return True
  return False
  
def assign_cause(row, dn_trans_vs_novel_df, transcripts_vs_assembly_dict, prot_vs_sample_proteins_dict, variant_annotation_dict, min_frac):
  """
  Given a row of the discrepancies
  table, assign a cause
  """
  cause = ''
  ## DN+/MTP- discrepancies
  if row['type'] == 1:
    # for unmatched genes
    if pd.isnull(row.iloc[2]):
      # check if transcript in novel sequences
      if query_in_psl(dn_trans_vs_novel_df, row.iloc[1], min_frac).empty:
        cause = 'gene_in_non-novel_region'
      else:
        cause = 'under_annotation'
    # for matched genes
    else:
      cause = 'read_mapping_issue'
  ## DN-/MTP+ discrepancies
  elif row['type'] == -1:
    # check if transcript is in assembly
    transcripts_vs_assembly_df = transcripts_vs_assembly_dict[row['sample']]
    # if transcript not found in assembly
    if query_in_psl(transcripts_vs_assembly_df, row.iloc[2], min_frac).empty:
      cause = 'incomplete_assembly'
    # if transcript found in assembly
    else:
      # check if protein is in sample
      prot_vs_sample_proteins_df = prot_vs_sample_proteins_dict[row['sample']]
      psl_hit = query_in_psl(prot_vs_sample_proteins_df, row.iloc[2], min_frac)
      # if protein is found
      if not psl_hit.empty:
        cause = 'clustering_issues'
      # if protein not found
      else:
        # check if evidence for pseudogenization exist
        sample_variant_annotation_dict = variant_annotation_dict[row['sample']]
        # if pseudogene - assign pseudogene cause
        if is_pseudogene(row.iloc[2], sample_variant_annotation_dict):
          cause = 'pseudogene'
        # if not, assign under-annotation cause
        else:
          cause = 'under_annotation'
  return cause


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Assign discrepancies causes')
  parser.add_argument('discrepancies_tsv', help='Path to discrepancies TSV')
  parser.add_argument('analysis_dir', help='Path to directory containing pipeline results')
  parser.add_argument('out_tsv', help='output discrepancies with assigned causes')
  parser.add_argument('--min_frac', type=float, default=0.9, help='Min fraction of query required to identify presence')

  args = parser.parse_args()

  # Read discrepancies table
  discrep_df = pd.read_csv(args.discrepancies_tsv, sep='\t')

  # Read psl file from DN_trans_vs_novel dir
  dn_trans_vs_novel_psl = os.path.join(args.analysis_dir, 'DN_trans_vs_novel/DN_trans_vs_novel.psl')
  dn_trans_vs_novel_df = parse_psl(dn_trans_vs_novel_psl)

  # Read psl files from transcripts_vs_assembly dir
  transcripts_vs_assembly_dir = os.path.join(args.analysis_dir, 'MTP_trans_vs_assembly')
  transcripts_vs_assembly_dict = create_data_dict(transcripts_vs_assembly_dir, 'psl', parse_psl)

  # Read psl files from prot_vs_sample_proteins dir
  prot_vs_sample_proteins_dir = os.path.join(args.analysis_dir, 'MTP_prot_vs_sample_proteins')
  prot_vs_sample_proteins_dict = create_data_dict(prot_vs_sample_proteins_dir, 'psl', parse_psl)

  # Read SnpEff BCFs
  variant_annotation_dir = os.path.join(args.analysis_dir, 'variant_calling')
  variant_annotation_dict =  create_data_dict(variant_annotation_dir, 'ann.bcf', parse_snpEff_bcf)

  # Assign causes
  apply_args = [dn_trans_vs_novel_df, transcripts_vs_assembly_dict, prot_vs_sample_proteins_dict, variant_annotation_dict, args.min_frac]
  discrep_df['cause'] = discrep_df.apply(assign_cause, args=apply_args, axis=1)

  # print output
  discrep_df.to_csv(args.out_tsv, sep='\t', index=False)
