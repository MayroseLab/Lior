from Bio import SeqIO
import numpy
import sys

def get_stats(fasta_path):
    results = {'Total length': 0,             
               'Total scaffolds': 0,
               '# of gaps': 0,
               '% gaps': 0,
               'N50': 0,
               'L50': 0,
               'N90': 0,
               'L90': 0,
               'Min scaffold length': 0,
               'Max scaffold length' :0
    }
    scaff_lengths = []
    with open(fasta_path) as f:
        for rec in SeqIO.parse(f,"fasta"):
            rec_length = len(rec)
            results['Total scaffolds'] += 1
            results['Total length'] += rec_length
            results['# of gaps'] += str(rec).count('N') + str(rec).count('n')
            scaff_lengths.append(rec_length)
    results['% gaps'] = results['# of gaps'] / results['Total length'] * 100
    results['Min scaffold length'] = min(scaff_lengths)
    results['Max scaffold length'] = max(scaff_lengths)
    sorted_scaff_lengths = sorted(scaff_lengths, reverse=True)
    csum = numpy.cumsum(sorted_scaff_lengths)
    csum2 = min(csum[csum >= sum(scaff_lengths)*0.5])
    results['L50'] = numpy.where(csum == csum2)[0][0]
    results['N50'] = sorted_scaff_lengths[results['L50']]
    csum9 = min(csum[csum >= sum(scaff_lengths)*0.9])
    results['L90'] = numpy.where(csum == csum9)[0][0]
    results['N90'] = sorted_scaff_lengths[results['L90']]
    return results

if __name__ == "__main__":
  res = get_stats(sys.argv[1])
  for s in ['Total length', 'Total scaffolds', '# of gaps', '% gaps', 'N50', 'L50', 'N90', 'L90', 'Min scaffold length', 'Max scaffold length']:
    print("{}\t{:,}".format(s,res[s]))

