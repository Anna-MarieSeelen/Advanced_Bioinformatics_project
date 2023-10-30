#!/usr/bin/env python3
"""
Author: Robert van der Klis
Student number: 1003241
Script to derive encoding used for fastq files
"""

# import statements
import gzip

# functions
def read_first_lines(n):
    """Returns the first n quality lines contained in a gzipped fastq file

    n: int, number of quality lines to read

    returns: list of strings, the first n quality lines
    """
    qualitylines = []
    for file in ['20180605.A-15R1_R1.fastq.gz', '20180613.A-26R5_R1.fastq.gz',
    '20180613.A-39R2_R1.fastq.gz', '20180605.A-33R5_R1.fastq.gz']:
        with gzip.open(f'../transcriptome/{file}', 'rb') as fo:
            for line in range(n):
                fo.readline()
                fo.readline()
                fo.readline()
                qualitylines.append(fo.readline().decode('utf-8').strip())
    return qualitylines

def count_qualities(qualitylines):
    """Counts the number of occurrences of each character in the quality lines

    qualitylines: list of strings, quality lines

    returns: dict as {character: number of occurrences}
    """
    qualitydict = {}
    for line in qualitylines:
        for ch in line:
            qualitydict.setdefault(ch, 0)
            qualitydict[ch] += 1
    return qualitydict

def main():
    """The main function of this module"""
    # Step 1: read quality lines
    qualitylines = read_first_lines(10000)
    
    # Step 2: count the number of occurrences of each quality
    qualdict = count_qualities(qualitylines)

    # Step 3: print to stdout
    print(qualdict)

if __name__ == "__main__":
    main()
