"""

"""

import networkx as nx
from networkx.algorithms import bipartite, matching
import sqlite3 as sql
import os
import sys
import pandas as pd
import numpy as np
from ete3 import Tree
from itertools import product, chain
from collections import Counter

def tidy_split(df, column, sep='|', keep=False):
  """
  Split the values of a column and expand so the new DataFrame has one split
  value per row. Filters rows where the column is missing.

  Params
  ------
  df : pandas.DataFrame
      dataframe with the column to split and expand
  column : str
      the column to split and expand
  sep : str
      the string used to split the column's values
  keep : bool
      whether to retain the presplit value as it's own row

  Returns
  -------
  pandas.DataFrame
      Returns a dataframe with the same columns as `df`.
  """
  indexes = list()
  new_values = list()
  df = df.dropna(subset=[column])
  for i, presplit in enumerate(df[column].astype(str)):
    values = presplit.split(sep)
    if keep and len(values) > 1:
      indexes.append(i)
      new_values.append(presplit)
    for value in values:
      indexes.append(i)
      new_values.append(value)
  new_df = df.iloc[indexes, :].copy()
  new_df[column] = new_values
  return new_df

def create_orthology_db(orthofinder_dir, weight_field='bitscore', overwrite=False):
  """
  Parse various OrthoFinder2 outputs to craeate
  an SQLite3 DB table with ortholog pairs. Each
  pair will have a score based on the field
  defined by weight_field (usually bitscore or
  pident), taken from Blast output as the 
  bi-directional mean (prot1->prot2 + prott2->prot1)/2.
  """
  # create DB
  db_path = os.path.join(orthofinder_dir, "orthologs.sqlite3")
  if os.path.exists(db_path):
    if not overwrite:
      print("Using existing DB %s" % db_path)
      return db_path
    elif overwrite:
      os.remove(db_path)
  con = sql.connect(db_path)
  cur = con.cursor()

  # Load Blast data
  #cur.execute("CREATE TABLE blast (qseqid text, sseqid text, pident real, length	integer, mismatch integer, gapopen	integer, qstart integer, qend integer, sstart integer, send integer, evalue real, bitscore real)")
  blast_dir = os.path.join(orthofinder_dir, 'WorkingDirectory')
  blast_files = [os.path.join(blast_dir, f) for f in os.listdir(blast_dir) if f.startswith('Blast')]
  colnames = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
  for bf in blast_files:
    df = pd.read_csv(bf, sep='\t', compression='gzip', header=None, names=colnames)
    df['qgenomeid'] = df['qseqid'].str.split('_', expand=True)[0]
    df['sgenomeid'] = df['sseqid'].str.split('_', expand=True)[0]
    df.to_sql('blast', con, if_exists='append', index=False)

  # Add protein names (by joining on protein IDs)
  cur.execute("CREATE TABLE seqnames (seqid text, seqname text)")
  seq_ids_file = os.path.join(blast_dir, 'SequenceIDs.txt')
  colnames = ['seqid','seqname']
  df = pd.read_csv(seq_ids_file, sep=': ', header=None, names=colnames)
  df.to_sql('seqnames', con, if_exists='append', index=False)
  cur.execute("CREATE TABLE blast_seqnames as SELECT b.*, s.* FROM (SELECT b.*, s.* FROM blast b INNER JOIN (SELECT seqid AS qseqid, seqname AS qseqname FROM seqnames) s ON b.qseqid = s.qseqid) b INNER JOIN (SELECT seqid AS sseqid, seqname AS sseqname FROM seqnames) s ON b.sseqid = s.sseqid")

  # Add genome names
  genome_ids_file = os.path.join(blast_dir, 'SpeciesIDs.txt')
  colnames = ['genomeid','genomename']
  df = pd.read_csv(genome_ids_file, sep=': ', header=None, names=colnames)
  df.to_sql('genomenames', con, if_exists='append', index=False)
  cur.execute("CREATE TABLE blast_seqnames_genomenames as SELECT b.*, s.* FROM (SELECT b.*, s.* FROM blast_seqnames b INNER JOIN (SELECT genomeid AS qgenomeid, genomename AS qgenomename FROM genomenames) s ON b.qgenomeid = s.qgenomeid) b INNER JOIN (SELECT genomeid AS sgenomeid, genomename AS sgenomename FROM genomenames) s ON b.sgenomeid = s.sgenomeid")

  # Load orthologs data
  orthologues_dir = [os.path.join(orthofinder_dir,d,'Orthologues') for d in os.listdir(orthofinder_dir)if d.startswith('Orthologues_')][0]
  orthologues_files = [os.path.join(dp, f) for dp, dn, fn in os.walk(orthologues_dir) for f in fn if 'Orthologues_' in dp and 'Xenologues' not in dp]
  for f in orthologues_files:
    df = pd.read_csv(f, sep='\t')
    df_t = df.iloc[:,[0,1,2]]
    df_t['prot1_genome'] = df_t.columns[1]
    df_t['prot2_genome'] = df_t.columns[2]
    df_t.columns = ['orthogroup','prot1', 'prot2', 'prot1_genome', 'prot2_genome']
    df_t = tidy_split(df_t, 'prot1', sep=', ')
    df_t = tidy_split(df_t, 'prot2', sep=', ')
    df_t['prot1'] = df_t['prot1'].str.replace('QI_','QI:').str.replace('AED_','AED:')
    df_t['prot2'] = df_t['prot2'].str.replace('QI_','QI:').str.replace('AED_','AED:')
    df_t.to_sql('orthologs', con, if_exists='append', index=False)

  # Take blast data for orthologs only
  cur.execute("CREATE TABLE blast_orthologs AS SELECT b.*, o.* FROM blast_seqnames_genomenames b INNER JOIN orthologs o ON b.qgenomename LIKE o.prot1_genome || '%' AND b.qseqname = o.prot1 AND b.sgenomename LIKE o.prot2_genome || '%' AND b.sseqname = o.prot2")
  # if there are multiple hits for the same proteins pair, choose the one with max weight_field
  cur.execute("CREATE TABLE blast_orthologs_max_weight AS SELECT *, max(%s) AS max_weight FROM blast_orthologs GROUP BY qseqid, sseqid" % weight_field)
  # calculate mean % identity of prot1->prot2 and prot2->prot1 and create final table
  cur.execute("CREATE TABLE blast_orthologs_bidirect AS SELECT orthogroup, prot1, prot2, prot1_genome, prot2_genome, avg(max_weight) OVER (PARTITION BY min(prot1, prot2), max(prot1, prot2), min(prot1_genome,prot2_genome), max(prot1_genome,prot2_genome)) weight FROM blast_orthologs_max_weight;")

  # finalize
  cur.close()
  con.close()
  return db_path
    
def tree_to_distance_matrix(tree, sister_only=False):
  """
  Takes an ete3 Tree object and returns
  a pandas DF with all pairwise distances
  between leaves.
  If sister_only, give distance inf to all
  non-sister leaves.
  """
  items = []
  genomes = []
  for g1 in tree:
    genome_name = g1.name
    d = []
    for g2 in tree:
      if sister_only and tree.get_distance(g1,g2,topology_only=True) > 1:
        dist = np.inf
      else:
        dist = tree.get_distance(g1,g2)
      d.append(dist)
    items.append((genome_name,d))
    genomes.append(genome_name)
  return pd.DataFrame.from_items(items, orient="index", columns=genomes)

def row_col_min(df, cutoff=0):
  """
  Find the row and column names of the
  min value in a DF. Ignore values < cutoff.
  Returns a (row,col) pair
  """
  # find min value
  vals = df.values
  vals[vals <= cutoff] = np.nan
  min_val = np.nanmin(vals)
  # find first row-column combination with min val
  for col in df:
    i = df.index[df[col] == min_val]
    if not i.empty:
      return (i[0],col)
  return (None,None)

def fill_list_to_length(x,l,v):
  """
  Add the value v to list x until
  it reaches length l
  """
  add_n = l - len(x)
  if add_n > 0:
    add = [v] * add_n
    x.extend(add)
  return x

class HomologyCluster(nx.Graph):

  """
  The class represents a raw orthogroup cluster
  from OrthoFinder2 Orthogroups.csv.
  In fact this is a homology group with paralogs
  present
  """

  def __init__(self, genomes, og_line, gene_tree=None):
    """
    Create an edge-less graph object, containing
    cluster genes as nodes with protin name as
    node name and genome in the 'genome' attribute.
    genomes - list of genome names in the
              order corresponding to the
              input line
    og_line - a line from the Orthogroups.csv
              file
    gene_tree - path to newick file with single
              gene tree of group genes
    """
    nx.Graph.__init__(self)
    self.has_paralogs = False
    fields = og_line.strip().split('\t')
    self.orthogroup = fields[0]
    for genome, genes in zip(genomes, fields[1:]):
      if genes == '':
        continue
      genes_list = genes.split(', ')
      if len(genes_list) > 1:
        self.has_paralogs = True
      for gene in genes_list:
        full_name = "%s_%s" %(genome, gene)
        full_name = full_name.replace('QI_','QI:').replace('AED_','AED:')
        self.add_node(full_name, genome=genome, gene_name=gene)
    self.gene_tree = None
    if gene_tree:
      self.gene_tree = Tree(gene_tree, format=1)
      for l in self.gene_tree:
        l.name = l.name.replace('QI_','QI:').replace('AED_','AED:').replace('_maker_proteins','.maker.proteins')

  def add_edges(self, db_con):
    """
    Add edges between cluster nodes, representing
    orthology inference, as given in the
    Orthologues files. Weights are assigned
    according to Blast % identity (mean of
    the two directions).
    All orthology data should be stored in an
    SQLite3 DB.
    """
    table_name = "blast_orthologs_bidirect"
    df = pd.read_sql_query("SELECT * FROM %s WHERE orthogroup = '%s'" %(table_name, self.orthogroup), con)
    df['prot1_full'] = df['prot1_genome'] + '_' + df['prot1']
    df['prot2_full'] = df['prot2_genome'] + '_' + df['prot2']
    df.apply(lambda g: self.add_edge(g['prot1_full'], g['prot2_full'], weight=g['weight']), axis=1)
    df.to_csv('tmp.tsv', sep='\t')

  def break_mwop(self, tree, allow_gene_copies=False):
    """
    Break a homology group into orthogroups
    where all proteins are orthologous to
    all proteins, using the
    Minimum Weight Orthogonal Partition (MWOP)
    criterion. A phylogenetic *species* tree
    is used to guide the order of actions in the
    algorithm. Optionally, the *gene* tree is used
    for keeping together recent gene copies.
    Returns a list of sets of gene names. Each
    set represents an orthogroup.
    See DOI:10.1007/978-3-642-23038-7_30 for
    details.
    """
    eps = 0.001
    tree = tree.copy()
    # find max number of genes in one genome
    #counts = [self.nodes[g]['genome'] for g in self.nodes]
    #counts = Counter(counts)
    #nm = counts.most_common(1)[0][1]

    # initialize gene sets
    genomes = [self.nodes[g]['genome'] for g in self.nodes]
    V = {genome: [] for genome in genomes}
    for gene in self.nodes:
      genome = self.nodes[gene]['genome']
      V[genome].append(set([gene]))
    # if allow_gene_copies mode is on, use gene tree
    # to cluster recent gene copies together
    if allow_gene_copies and self.gene_tree:
      # create a list of all sister genes
      sisters = set()
      for g in self.gene_tree:  # g is a leaf (gene)
        sis = g.get_sisters()
        for s in sis:
          if s.is_leaf() and (s.name,g.name) not in sisters:
            sisters.add((g.name,s.name))
      # only keep sister pairs from the same genome
      gene_copies = set()
      for sis_pair in sisters:
        if self.nodes[sis_pair[0]]['genome'] == self.nodes[sis_pair[1]]['genome']:
          gene_copies.add(sis_pair)
      # create sets of gene copies and remove single-gene sets
      rm = set(chain.from_iterable(gene_copies))
      for genome in V:
        new_Vi = []
        for s in V[genome]:
          if s == set() or next(iter(s)) not in rm:
            new_Vi.append(s)
        V[genome] = new_Vi
      for cp in gene_copies:
        cp_genome = self.nodes[cp[0]]['genome']
        V[cp_genome].append(set(cp))
    # complete with empty sets
    nm = max([len(l) for l in V.values()])
    for genome in V:
      V[genome] = fill_list_to_length(V[genome], nm, set())

    # remove tree leaves with no genes
    for genome in tree:
      if genome.name not in V:
        genome.delete()

    # traverse tree and create new orthogroups
    n_leaves = len(tree.get_tree_root())
    while n_leaves > 1:
      # find two closest genomes based on tree
      dist_matrix = tree_to_distance_matrix(tree, sister_only=True)
      closest_genomes = row_col_min(dist_matrix)
      # create bipartite graph (BG)
      # (BG has integers as labels, so Vi_d and Vj_d
      # store the integer --> gene set mapping)
      Vi = V[closest_genomes[0]]
      Vi_d = {n: Vi[n] for n in range(nm)}
      Vj = V[closest_genomes[1]]
      Vj_d = {n: Vj[n-nm] for n in range(nm,nm*2)}
      bipartite_graph = bipartite.complete_bipartite_graph(nm,nm)
      # assign weights to BG edges
      for edge in bipartite_graph.edges:
        # look for corresponding edge in homology graph (HG)
        # first, translate BG integers to gene sets
        g1 = Vi_d[edge[0]]
        g2 = Vj_d[edge[1]]
        # calculate edge weights as mean of weights between sets
        pairs = product(g1,g2)
        total_weight = 0
        pairs_involved = 0
        for p in pairs:
          if p in self.edges: # orthologs
            total_weight += self.edges[p]['weight']
          else: # paralogs
            total_weight += -eps
          pairs_involved += 1
        if pairs_involved == 0: # when using an empty set
          bipartite_graph.edges[edge]['weight'] = eps
        else:
          bipartite_graph.edges[edge]['weight'] = total_weight/pairs_involved
      # find BG maximum weight matching (MWM)
      # ensure matches are always from Vi to Vj
      mwm = matching.max_weight_matching(bipartite_graph,maxcardinality=True)
      mwm = {(min(x),max(x)) for x in mwm}
      # update gene sets according to MWM
      # this is done by set unions on every matching
      new_Vi = []
      for match in mwm:
        g1 = Vi_d[match[0]]
        g2 = Vj_d[match[1]]
        g1 = g1.union(g2)
        new_Vi.append(g1)
      V[closest_genomes[0]] = new_Vi
      # then remove other genome from V
      del(V[closest_genomes[1]])
      # update tree - remove leaf corresponding to Vj
      # but update branch length of Vi to mean of branch
      # lengths of Vi and Vj
      Vj_bl = (tree&closest_genomes[1]).dist
      j = tree.search_nodes(name=closest_genomes[1])[0]
      j.delete()
      Vi_bl = (tree&closest_genomes[0]).dist
      new_Vi_bl = (Vi_bl + Vj_bl)/2
      (tree&closest_genomes[0]).dist = new_Vi_bl
      # calculate number of leaves left
      n_leaves = len(tree.get_tree_root())

    # finish when all leaves were merged
    return V

if __name__ == "__main__":
  orthofinder_dir = sys.argv[1]
  weight_field = sys.argv[2] # e.g. pident, bitscore
  orthogroups_out = sys.argv[3]
  
  # create orthologs DB
  db_path = create_orthology_db(orthofinder_dir)
  con = sql.connect(db_path)

  # species tree
  orthologues_dir = [os.path.join(orthofinder_dir,d) for d in os.listdir(orthofinder_dir)if d.startswith('Orthologues_')][0]
  tree_file = os.path.join(orthologues_dir, 'SpeciesTree_rooted_node_labels.txt')
  tree = Tree(tree_file, format=1)

  # Calculate orthogroups
  orthogroups_file = os.path.join(orthofinder_dir, 'Orthogroups.csv')
  gene_trees_dir = [os.path.join(orthofinder_dir,d,"Recon_Gene_Trees") for d in os.listdir(orthofinder_dir) if d.startswith('Orthologues_')][0]

  with open(orthogroups_file) as f, open(orthogroups_out, 'w') as fo:
    genomes = f.readline().strip().split('\t')
    for line in f:
      og = line.split('\t')[0]
      print(og)
      og_tree = os.path.join(gene_trees_dir,"%s_tree.txt" % og)
      if os.path.isfile(og_tree):
        hc = HomologyCluster(genomes, line, gene_tree=og_tree)
      else:
        hc = HomologyCluster(genomes, line)
      if hc.has_paralogs and len(hc.nodes) > 3:
        hc.add_edges(con)
        orthogroups = hc.break_mwop(tree, allow_gene_copies=True)
        orthogroups = list(orthogroups.values())[0]
      else:
        orthogroups = [list(hc.nodes)]
      #if og == "OG0002703":
      #  print(orthogroups)
      for og_break in orthogroups:
        if len(og_break) > 0:
          print('\t'.join(og_break), file=fo)
