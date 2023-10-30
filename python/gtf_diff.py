"""
Finds gene-ID's that are only expressed in the reference-mapping and not in
the query-mapping. By giving -g option, the script can compare gene-ID's
in two reference genomes. The genes are written to a file, line by line,
and the amount can be assesed by using wc.

Author: JN Bondevik

Usage: python3 gtf_diff.py [-g] ref_gtf query_gtf output
-g [optional]: for comparing two unmapped reference genomes
"""

from sys import argv
import os
import re

def check_args():
    """Checks if required command line arguments are given. If -g is not given
    it is set as False in the argv-list."""
    if len(argv) < 4:
        raise ValueError('Not given enough arguments.\n'\
        'Usage: python3 gtf_diff.py [-g] ref_gtf query_gtf output.\n'\
        '-g: for comparing two unmapped reference-genomes.\n')
    elif len(argv) < 5:
        argv.insert(1, None)
    elif len(argv) == 5:
        if argv[1] != '-g':
            raise ValueError(f'"{argv[1]}" is not recognized as an argument.')
    return argv[1], argv[2], argv[3], argv[4]

def list_samples(path):
    """Finds and adds all unique gtf-filepaths in a directory to a set."""
    samples = {path+sample for sample in os.listdir(path) if "gtf" in sample}
    return samples

def parse_genes(gtf_file_path, g):
    """Parse one gtf-file and return all unique gene-ID's.

    gtf (str): path to a gtf-file.
    g (bool): an option to specify if coverage should be checked.

    Returns a set with all unique gene-ID's with coverage from the gtf.
    If -g is True, the coverage is not considered."""
    genes = set()
    with open(gtf_file_path, 'r') as gtf:
        for line in gtf:
            if not line.startswith('#'):
                if g:
                    name = re.match(r'.*gene_id "([^"]+)"', line).group(1)
                    if name:
                        genes.add(name)
                else:
                    match = re.match(r'.*gene_id "([^"]+)".*cov "([^"]+)', line)
                    if match:
                        name = match.group(1)
                        cov = float(match.group(2))
                        if cov > 0:
                            genes.add(name)
    return genes

def merge_genes(gtf_path, g):
    """Merges all unique gene-ID's from all gtf-files found in a directory.
    The gene-ID's are added to a set.

    gtf_path: path to the directory with the gtf-files.

    Returns a set of all genes found in all gtf-files in a directory."""
    genes = set()
    for filepath in list_samples(gtf_path):
        for gene in parse_genes(filepath, g):
            genes.add(gene)
    return genes

def compare_genes(ref_genes, query_genes):
    """Compares reference-genes to query-genes. If a reference gene is not 
    found in the query, it is added to a set of unique genes.

    ref_genes (set): genes found in a reference gtf-file or -directory.
    query_genes (set): genes found in a query gtf-file or -directory. 

    Returns a set of reference genes not found in the query genes."""
    unique = {gene for gene in ref_genes if not gene in query_genes}
    return unique

def write_unique(unique_genes, output):
    """Writes all genes in unique genes to an output.
    
    unique_genes (set): all genes to write to output.
    output (string): name of the output-file to create.

    Returns the name of the output-file."""
    with open(output, 'w') as out:
        for gene in unique_genes:
            out.write(f"{gene}\n")
    return output

def main():
    # Step 1: Assign command line arguments
    g_opt, ref, query, output  = check_args()
    # Step 2: Check if path is to a file or not
    if os.path.isfile(ref):
        # Step 3a: Compare the gene-ID's
        unique_genes = compare_genes(parse_genes(ref, g_opt), parse_genes(query, g_opt))
    else:
        # Step 3b: Compare the merged gene-ID's
        unique_genes = compare_genes(merge_genes(ref, g_opt), merge_genes(query, g_opt))
    # Step 4: Write to output
    write_unique(unique_genes, output)
    print("Complete. ")

if __name__ == "__main__":
    main()
