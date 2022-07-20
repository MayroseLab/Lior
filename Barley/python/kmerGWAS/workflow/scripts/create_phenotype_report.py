import sys
import os
import re
from intervaltree import IntervalTree, Interval
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import plotly.io as pio
pio.templates.default = "plotly_white"

pheno_name = sys.argv[1]
pheno_tsv = sys.argv[2]
signif_kmers_by_acc_tsv = sys.argv[3]
kmer_mapping_tsv = sys.argv[4]
ref_name = sys.argv[5]
chr_len_tsv = sys.argv[6]
ref_gff = sys.argv[7]
contig_mapping_dir = sys.argv[8]
genomes_tsv = sys.argv[9]
min_chr_len = int(sys.argv[10])
n_top_acc = int(sys.argv[11])
out_html = sys.argv[12]

### FUNCTIONS

def manhattan_plot():
    fig = px.scatter(mh_plot_df, x='Pos', y='-logP', color='map_type',
                 facet_col='Chr', facet_col_spacing=0.01,
                 width=2000, height=1000,
                labels={'Pos':''},
                    title='k-mer GWAS Manhattan plot')
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1], y=-0.12))

    fig.update_xaxes(showgrid=False, showticklabels=False, showline=True, mirror=True, linecolor='black',
                     zeroline=False, ticks="", matches=None)
    fig.update_yaxes(showgrid=False, showline=True, mirror=True, linecolor='black', matches=None)
    fig.update_layout(legend_title_text='')
    chr_list_cp = chr_list[:]
    def update_xaxis(xaxis):
        chrom = chr_list_cp.pop(0)
        xaxis.update(range=[1, chr_len_df.loc[chr_len_df['Chr'] == chrom]['len']], title=chrom)
    fig.for_each_xaxis(update_xaxis)
    return fig

def fix_attr(attr_list):
    fixed = []
    for a in attr_list:
        if '=' in a:
            fixed.append(a)
        else:
            fixed[-1] += ';%s' % a
    return fixed

def gff_to_df(gff):
    gff_headers = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff_df = pd.read_csv(gff, sep='\t', comment='#', names=gff_headers)
    gff_df = gff_df.query('type == "gene"')
    gff_df['attributes_dict'] = gff_df['attributes'].apply(lambda s: {a.split('=')[0]: a.split('=')[1] for a in fix_attr(s.split(';'))})
    gff_df['gene_ID'] = gff_df['attributes_dict'].apply(lambda x: x['ID'])
    gff_df = gff_df[['seqid', 'start', 'end', 'gene_ID',]]
    return gff_df

def gene_trace(gene_row):
    x = [gene_row['start'], gene_row['end']]
    y = [0,0]
    return go.Scatter(x=x, y=y, mode='lines', line={'color': 'blue', 'width': 4},
                     hovertext=gene_row['gene_ID'])

def list_to_pairs(l):
    it = iter(l)
    return list(zip(it,it))

def cigar_to_ivtree(cigar_str):
    cigar_pairs = list_to_pairs(re.split(r'([SHMID])', cigar_str))
    cigar_ivt = IntervalTree()
    s = 1
    for p in cigar_pairs:
        if p[1] not in 'SHMID':
            continue
        l = int(p[0])
        cigar_ivt[s:s+l] = p[1]
        s += l
    return cigar_ivt

def cigar_trace(cigar_ivt, start_from=1, name=''):
    x = []
    y = []
    c = start_from
    for iv in sorted(cigar_ivt):
        iv_len = iv.end - iv.begin
        if iv.data == 'M':
            x.append(c)
            y.append(0)
            x.append(c + iv_len)
            y.append(0)
        elif iv.data in {'S','H'}:
            x.append(c)
            y.append(1)
            x.append(c + iv_len)
            y.append(1)
        elif iv.data == 'I':
            x.append(c - iv_len/2)
            y.append(1)
            x.append(c + iv_len/2)
            y.append(1)
        elif iv.data == 'D':
            x.append(c)
            y.append(-1)
            x.append(c + iv_len)
            y.append(-1)
        if iv.data in {'M','S','H','D'}:
            c += iv_len
    return go.Scatter(x=x, y=y, mode='lines', hovertext=name)

def sam_to_df(sam_p, header_pref=''):
    sam_headers = ["QNAM", "FLAG", "RNAM", "POS", "MAPQ", "CIGAR", "RNEX", "PNEX", "TLEN", "SEQ", "QUAL"]
    df = pd.read_csv(sam_p, sep='\t', comment='@', usecols=range(11), names=sam_headers)
    df = df[["QNAM", "RNAM", "POS", "CIGAR"]]
    df.columns = ["%s_%s" %(header_pref, col) for col in df.columns]
    return df

def create_chrom_plot(chrom=''):
    # Create plot with subplots
    n_proj_acc = len(topp_acc_list)
    n_subplots = 2+n_proj_acc*2
    fig = make_subplots(rows=n_subplots, cols=1, shared_xaxes=True,
                   row_heights=[100,50]*n_proj_acc + [100,20],
                   row_titles=sub_plot_titles)
           
    # Reference genes
    chrom_genes_df = ref_gff_df.query('seqid == @chrom')
    for gene in chrom_genes_df.iterrows():
        fig.add_trace(gene_trace(gene[1]), row=n_subplots, col=1)
    
    # Direct mapping
    chrom_direct_map_df = kmer_mapping['Direct'].query('Chr == @chrom')
    fig.add_trace(go.Scatter(x=chrom_direct_map_df['Pos'], y=chrom_direct_map_df['-logP'], mode='markers', marker={'size':2}),
              row=n_subplots-1, col=1)
    
    # k-mer projections
    i = 1
    for acc_name in topp_acc_list:
        chrom_kmers_proj_df = kmer_mapping[acc_name].query('Chr == @chrom')
        fig.add_trace(go.Scatter(x=chrom_kmers_proj_df['Pos'], y=chrom_kmers_proj_df['-logP'], mode='markers', marker={'size':2}),
                        row=i, col=1)
        i += 2
    
    # contig mapping
    i = 2
    for acc_name in topp_acc_list:
        chrom_contigs_map_df = contig_mapping[acc_name].query('contig_RNAM == @chrom')
        for row in chrom_contigs_map_df.groupby(['contig_QNAM','contig_POS']).first().reset_index().iterrows():
            row = row[1]
            fig.add_trace(cigar_trace(row['contig_CIGAR_ivt'], start_from=row['contig_POS'], name=row['contig_RNAM']), row=i, col=1)
        i += 2
     
    # Layout
    total_height = 350 + 150*n_proj_acc
    fig.update_layout(showlegend=False, width=2000, height=total_height, title='Chromosome %s' % chrom, margin={'l': 200})
    for yaxis_n in range(2,n_subplots+1,2):
        fig.layout['yaxis%s' % yaxis_n]['range'] = [-2,2]
        fig.layout['yaxis%s' % yaxis_n]['showticklabels'] = False
        fig.layout['yaxis%s' % yaxis_n]['ticks'] = ''
        fig.layout['yaxis%s' % yaxis_n]['fixedrange'] = True
    fig.update_yaxes(mirror=True)
    fig.update_xaxes(mirror=True)
    for annotation in fig['layout']['annotations']: 
        annotation['textangle']=0
        annotation['x'] = -0.06
        annotation['font'] = {'size': 12}

    return fig

def write_figs_to_html(out_html, figs_list):
    if os.path.exists(out_html):
        os.remove(out_html)
    with open(out_html, 'a') as fo:
        for fig in figures:
            print(fig.to_html(full_html=False, include_plotlyjs='cdn'), file=fo)

### MAIN

if __name__ == "__main__":

    figures = []

    # Phenotype histogram
    pheno_df = pd.read_csv(pheno_tsv, sep='\t')
    n_accessions = pheno_df.shape[0]
    plot_title = "%s (N=%s) - phenotype values" %(pheno_name.replace('_',' '), n_accessions)
    max_val = pheno_df['phenotype_value'].max()
    pheno_hist_fig = px.histogram(pheno_df['phenotype_value'], title=plot_title)
    pheno_hist_fig = pheno_hist_fig.update_xaxes(range=[0, max_val])
    pheno_hist_fig = pheno_hist_fig.update_layout(showlegend=False)
    figures.append(pheno_hist_fig)

    # Udi plot
    signif_kmers_by_acc_df = pd.read_csv(signif_kmers_by_acc_tsv, sep='\t')
    signif_kmers_by_acc_df.sort_values(by='kmers_count', inplace=True, ascending=False)
    signif_kmers_by_acc_df.columns = ['Accession', 'kmers_count', pheno_name]
    udi_plot_fig = px.bar(signif_kmers_by_acc_df, y='kmers_count', color=pheno_name,
                     title='Udi plot', hover_data=['Accession', 'kmers_count', pheno_name])
    udi_plot_fig = udi_plot_fig.update_xaxes(title='')
    udi_plot_fig = udi_plot_fig.update_yaxes(title='Significant k-mers')
    udi_plot_fig = udi_plot_fig.update_layout(legend_title_text=pheno_name)
    figures.append(udi_plot_fig)

    
    topp_acc_list = list(signif_kmers_by_acc_df.iloc[:n_top_acc,0])
    # remove top accessions if an assembly does not exist
    genomes_df = pd.read_csv(genomes_tsv, sep='\t')
    lq_genomes_with_assembly = set(genomes_df.query('genome_type == "LQ"')['genome_name'])
    topp_acc_list = [acc_name for acc_name in topp_acc_list if acc_name in lq_genomes_with_assembly]

    # Manhattan plot
    kmer_mapping_df = pd.read_csv(kmer_mapping_tsv, sep='\t', index_col=0)
    if kmer_mapping_df.empty:
        write_figs_to_html(out_html, figures)
        print("No significant k-mers. Terminating.")
        exit(0)

    chr_len_df = pd.read_csv(chr_len_tsv, sep='\t', names=['Chr','len'])
    chr_list = sorted(list(chr_len_df.query('len >= @min_chr_len')['Chr']))
    mh_direct_map_df = kmer_mapping_df.loc[~ kmer_mapping_df['Chr'].isna()][['Chr', 'Pos', '-logP']]
    mh_direct_map_df['map_type'] = 'Direct'
    if topp_acc_list:
        acc_top_kmers = topp_acc_list[0]
        mh_proj_map_df = kmer_mapping_df.loc[~ kmer_mapping_df['Chr_%s' % acc_top_kmers].isna()][['Chr_%s' % acc_top_kmers, 'Pos_%s' % acc_top_kmers, '-logP']]
        mh_proj_map_df.columns = ['Chr', 'Pos', '-logP']
    else:
        mh_proj_map_df = pd.DataFrame({'Chr':[], 'Pos':[], '-logP':[]})
    mh_proj_map_df['map_type'] = 'Projection'
    mh_plot_df = pd.concat([mh_direct_map_df, mh_proj_map_df])
    mh_plot_df = mh_plot_df.loc[(mh_plot_df['Chr'].isin(chr_list)) | (mh_plot_df['Chr'].isna())]
    mh_plot_df.sort_values(by='Chr', inplace=True)
    mh_plot_fig = manhattan_plot()
    figures.append(mh_plot_fig)

    # Chromosome figures
    # create subset DFs
    kmer_mapping = {}
    kmer_mapping['Direct'] = kmer_mapping_df.loc[~ kmer_mapping_df['Chr'].isna()][['Chr','Pos','-logP']]

    for acc_name in topp_acc_list:
        acc_df = kmer_mapping_df.loc[~ kmer_mapping_df['Chr_%s' % acc_name].isna()][['Chr_%s' % acc_name,'Pos_%s' % acc_name, '-logP']]
        acc_df.columns = ['Chr','Pos','-logP']
        kmer_mapping[acc_name] = acc_df

    # locate and parse relevant contig mapping files
    contig_mapping = {}
    for acc_name in topp_acc_list:
        contig_sam = os.path.join(contig_mapping_dir, 'pass_threshold_5per_%s.signif_contigs_vs_%s.filter.sam' %(acc_name, ref_name))
        contig_df = sam_to_df(contig_sam, header_pref='contig')
        contig_df['contig_CIGAR_ivt'] = contig_df['contig_CIGAR'].apply(cigar_to_ivtree)
        contig_mapping[acc_name] = contig_df

    ref_gff_df = gff_to_df(ref_gff)

    # Sub plot titles
    sub_plot_titles = []
    for acc_name in topp_acc_list:
        sub_plot_titles.append('%s<br>k-mers' % acc_name)
        sub_plot_titles.append('%s<br>contigs' % acc_name)
    sub_plot_titles.extend(['Direct<br>mapping', '%s<br>genes' % ref_name])

    # Create figure per chrom
    for chrom in chr_list:
        figures.append(create_chrom_plot(chrom=chrom))

    ### Write to HTML
    write_figs_to_html(out_html, figures)
