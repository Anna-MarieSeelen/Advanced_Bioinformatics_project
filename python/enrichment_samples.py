#!/usr/bin/env python3
"""
Author: Robert van der Klis
Student number: 1003241
Script to find samples in the dir containing DE output, and put them in a txt
Usage: python3 enrichment_samples output_dir
output_dir: string, name of the directory containing DE analysis output
"""

# import statements
import subprocess
from sys import argv

# functions
def list_files(dir):
    """Lists files in directory
    
    dir: string, directory to look in
    
    output: list of strings, contents of current dir"""
    contents = subprocess.check_output(f'ls ./{dir}', shell=True)
    contents = contents.decode('utf-8').split()
    return contents

def main():
    """The main function of this module"""
    # Step 1: return contents of DE output dir
    contents = list_files(argv[1])

    # Step 2: write contents of DE output dir to enrichment_samples.txt
    with open('sample_list.txt', 'w') as fo:
        [fo.write(argv[1] + '/' + sample + '\n') for sample in contents]

if __name__ == "__main__":
    main()