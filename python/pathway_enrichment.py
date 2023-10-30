#!/usr/bin/env/python3
"""
Author: Robert van der Klis
Student number: 1003241
Script to do pathway enrichment analysis on Arabidopsis Thaliana
Usage: python3 pathway_enrichment.py aracyc_pathways sample_list
aracyc_pathways: tsv file of aracyc arabidopsis pathways with genes
sample_list: txt file of all .tsv files with DEG counts to include in analysis
separated by newlines
"""

# import statements
from scipy import stats
from sys import argv

# functions
def sample_list(sample_file):
    """Reads file containing samples to analyse
    
    sample_file: string, name of the file containing the samples to analyse
    
    output: list of strings, the names of the samples
    """
    with open(sample_file) as fo:
        samples = fo.readlines()
    samples = [sample.strip() for sample in samples]
    return samples

def read_tsv(tsvfile):
    """Parse tab-separated file
    
    tsvfile: string, name of the file containing the pathways
    
    output: list, each element is one line in the input file, split on tabs
    """
    with open(tsvfile) as fo:
        lines = fo.readlines() 
    output = []
    for line in lines:
        words = line.split('\t')
        output.append([word.strip() for word in words])
    return output

def split_genes(diff_exp):
    """Split genes into > 0 log2fold, and < 0 log2fold
    
    diff_exp: list, containing the lines of the DE analysis, split on tabs
    
    output:
    upreg: list, containing the lines with > 0 log2fold
    downreg: list, containing the lines with < 0 log2fold
    """
    upreg, downreg = [], []
    for line in diff_exp[1:]:
        if line[2] == 'NA':
            continue
        if float(line[2]) > 0:
            upreg.append(line)
        else:
            downreg.append(line)
    return upreg, downreg

def parse_pathways(pathways):
    """Parses lines of the pathways file into a dict
    
    pathways: list of lines, containing arabidopsis pathways
    
    output: dict, as {gene: pathway}
    """
    pathwaydict = {}
    for row in pathways:
        # Filter out pathways with 'unknown' as gene id
        if row[-2] != 'unknown':
            pathwaydict[row[-2]] = row[1]
    return pathwaydict

def count_pathway_degs(diff_exp, pathwaydict):
    """Counts differentially expressed genes per pathway
    
    diff_exp: list of lines, differential expression analysis results
    pathwaydict: dict, as {gene: pathway}
    
    output: 
    pathway_counts: dict, as {pathway: [deg, non-deg]}. DEGs per pathway
    total_counts: list, as [deg, non-deg]. Total DEGs
    """
    pathway_counts = {}
    total_counts = [0, 0]
    for row in diff_exp[1:]:
        geneid = row[0].split('|')[0]
        padj = float(row[-1])
        if padj < 0.05:
            total_counts[0] += 1
        else:
            total_counts[1] += 1
        try:
            # Pathwaydict contains {gene: pathway}, so create pathway_counts 
            # by looking for each gene, using the pathway as
            # key, and creating 2 counters (one for deg, one for non-deg)
            pathway_counts.setdefault(pathwaydict[geneid], [0, 0])
            if padj < 0.05:
                pathway_counts[pathwaydict[geneid]][0] += 1
            else:
                pathway_counts[pathwaydict[geneid]][1] += 1
        # If some genes are not in the original pathways file, simply skip. They
        # are added to total_counts, however.
        except KeyError:
            continue
    return pathway_counts, total_counts

def fisher_test(pathway_counts, total_counts):
    """Tests each pathway for significance
    
    pathway_counts: dict, as {pathway: [deg, non-deg]}. DEGs per pathway
    total_counts: list, as [deg, non-deg]. Total DEGs
    
    output: dict, as {pathway: [odds ratio, p value]}
    """
    nonsig_pathways, sig_pathways = {}, {}
    for pathway, currcounts in pathway_counts.items():
        restcounts = [total - count for total, count in zip(total_counts, currcounts)]
        # do one-sided test
        pval = stats.fisher_exact([currcounts, restcounts], alternative='greater')[1]
        if pval > 0.05:
            nonsig_pathways[pathway] = pval
        else:
            sig_pathways[pathway] = pval
    return nonsig_pathways, sig_pathways

def main():
    """The main function of this module"""
    # Initialize dicts for pathways
    n_up_sig_pathways = {}
    n_up_sig_pathways.setdefault('araport', {})
    n_up_sig_pathways.setdefault('tair', {})
    n_down_sig_pathways = {}
    n_down_sig_pathways.setdefault('araport', {})
    n_down_sig_pathways.setdefault('tair', {})

    # Step 1: get the pathways, and get the list of sample files
    pathways = read_tsv(argv[1])
    samplist = sample_list(argv[2])

    # Step 2: put pathways into dict
    pathwaydict = parse_pathways(pathways)

    for sample in samplist:
        # Step 3: read the genes expression, split on tabs
        diff_exp = read_tsv(sample)
        
        # Step 4: split genes into log2fold > 0 and log2fold < 0
        upreg, downreg = split_genes(diff_exp)

        # Step 5: count DEGs, non-DEGs per pathway, and in total
        up_degs = count_pathway_degs(upreg, pathwaydict)
        down_degs = count_pathway_degs(downreg, pathwaydict)
        up_pathway_counts, up_total_counts = up_degs
        down_pathway_counts, down_total_counts = down_degs
        
        # Step 6: statistical analysis on DEGs, non-DEGs per pathway, vs total
        up_test = fisher_test(up_pathway_counts, up_total_counts)
        down_test = fisher_test(down_pathway_counts, down_total_counts)
        up_sig_pathways = up_test[1]
        down_sig_pathways = down_test[1]

        # Step 7: print results
        # Print up- and downregulated pathways per sample
        print('\n\n' + sample)
        print('up:')
        [print(key, end=', ') for key in up_sig_pathways]
        print()
        print('down:')
        [print(key, end=', ') for key in down_sig_pathways]

        # Code necessary to print the most ubiquitously enriched pathways
        for pathway in up_sig_pathways.keys():
            n_up_sig_pathways['araport'].setdefault(pathway, 0)            
            n_up_sig_pathways['tair'].setdefault(pathway, 0)
            if 'expressions_araport11' in sample:
                n_up_sig_pathways['araport'][pathway] += 1
            else:
                n_up_sig_pathways['tair'][pathway] += 1

        for pathway in down_sig_pathways.keys():
            n_down_sig_pathways['araport'].setdefault(pathway, 0)            
            n_down_sig_pathways['tair'].setdefault(pathway, 0)
            if 'expressions_araport11' in sample:
                n_down_sig_pathways['araport'][pathway] += 1
            else:
                n_down_sig_pathways['tair'][pathway] += 1

    # Step 8: print the pathways enriched in at least n samples
    n = 6
    print()
    print('up:')
    for annotation, pathways in n_up_sig_pathways.items():
        for pathway, num_samples in pathways.items():
            if num_samples >= n:
                print(annotation, pathway, num_samples)

    print('down:')
    for annotation, pathways in n_down_sig_pathways.items():
        for pathway, num_samples in pathways.items():
            if num_samples >= n:
                print(annotation, pathway, num_samples)

if __name__ == "__main__":
    main()