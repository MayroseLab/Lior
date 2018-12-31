from __future__ import print_function, division
import argparse
from gff3 import Gff3
import pandas as pd


### FUNCTIONS

def mrna_aed(gff):
    """
    Takes a gff3 object and returns a dictionary
    of AED scores with mRNA names as keys.
    """
    return ({mrna_id: float(gff.features[mrna_id][0]['attributes']['_AED']) for mrna_id in gff.features if gff.features[mrna_id] and gff.features[mrna_id][0]['type'] == "mRNA"}

def mrna_utr(gff):
    """
    Takes a gff3 object and returns a dictionary.
    Keys are mRNA IDs and values are UTR status
    -1 - illegal UTR; 0 - no UTRs; 1 - missing 3' or 5' UTR; 2 - legal UTRs
    """
    res = {}
    translate = {'five_prime_UTR':'5', 'three_prime_UTR':'3', 'CDS': 'C'}
    mrna_lines = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        mrna_id = mrna_line['attributes']['ID']
        descendants = gff.descendants(mrna_id.line_index)
        descendants_features_order = [f['type'] for f in descendants if f['type'] in ['five_prime_UTR', 'CDS', 'three_prime_UTR']]
        descendants_features_order_simp = ''.join([translate[t] for t in descendants_features_order])
        if '5' not in descendants_features_order_simp and '3' not in descendants_features_order_simp:
            stat = 0
        elif 'C5C' in descendants_features_order_simp or 'C3C' in descendants_features_order_simp:
            stat = -1
        elif '5' not in descendants_features_order_simp or '3' not in descendants_features_order_simp:
            stat = 1
        else:
            stat = 2
        res[mrna_id] = stat
    return res

def prot_busco(busco_full):
    """
    Reads a BUSCO full report and returns a dict
    with proteins in which BUSCOs were found. {gene name: BUSCO ID}
    """
    res = {}
    with open(busco_full):
        for line in f:
            fields = line.split('\t')
            if fields[1] == "Complete" or fields[1] == "Duplicated":
                res[fields[2]] = fields[0]
    return res

def prot_similarity(blast_res):
    """
    Reads a blast tsv result and returns a dict with
    gene names as keys and query coverage as values.
    Assumes query name at column 2 and query coverage at column 10.
    """
    res = {}
    with open(blast_res) as f:
        for line in f:
            fields = line.strip().split('\t')
            res[fields[1]] = fields[9]
    return res

### MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_gff', help="Input gff3 annotation file")
    parser.add_argument('--busco_result', default=None, help="BUSCO full report output")
    parser.add_argument('--ips_result', help="InterProScan output")
    parser.add_argument('--blast_result', help="Blast search result of proteins vs. DB")
    parser.add_argument('out_report', help="Output QA report")
    args = parser.parse_args()

    gff_obj = Gff3(args.in_gff)
    qa_methods = [("AED",mrna_aed,gff_obj), ("UTR",mrna_utr,gff_obj)]
    if args.busco_result:
        qa_methods.append(('BUSCO',prot_busco,args.busco_result))
    if args.blast_result:
        qa_methods.append(('BLAST',prot_similarity,args.blast_result))
    

    qa_data = [ pd.Series(m[1].__call__(m[2]), name=m[0]) for m in qa_methods ]
    qa_df = pd.concat(qa_data, axis = 1)
    with open(args.out_report,'w') as fo
        print(qa_df.to_string(), file=fo)
