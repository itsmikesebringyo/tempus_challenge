# Mike Sebring
# Tempus Coding Challenge
# October 19, 2020

import requests
import json
import os


filepath = 'Challenge_data_(1).vcf'


#### test with pandas read
import pandas as pd

with open(filepath) as ft:
    start = ft.read().find('#CHROM')
    ft.seek(start+1)
    df1 = pd.read_csv(ft, sep='\t')


df1.columns
df1
df1['CHROM'][2]

df1['variant_id'] = df1.agg(lambda x: f"{x['CHROM']}-{x['POS']}-{x['REF']}-{x['ALT']}", axis=1)
df1['variant_id'][1000]

###### helper functions to grab the necessary values

def get_allele_freq(variant_id):
    """Takes variant ID as input and returns allele frequency."""
    
    # Hit the ExAC API and convert the results into JSON.
    variant_info = requests.get(f'http://exac.hms.harvard.edu/rest/variant/variant/{variant_id}').json()
    
    # Try getting the allele frequency from the API result if possible.
    # Otherwise, we need to use the allele frequency provided in in source file.
    try:
        allele_frequency = variant_info['allele_freq']
    except:
        allele_frequency = 0.5

    return allele_frequency

test_allele_freq = get_allele_freq(df1['variant_id'][1000])
print(test_allele_freq)

def create_info_dict(info_col):
    """Takes INFO column of variant and returns dictionary of keys/values."""

    # First, we need to split the entire column so the key/value pairings are
    # separated from each other.
    info_split = info_col.split(';')

    # Next, we iterate through each key/value pairing and write them 
    # into our INFO dictionary.
    info_dict = {}
    for pairing in info_split:
        key_val = pairing.split('=')
        
        # Include a try/except clause so that we can evaluate int and
        # float variables but not error out on str. 
        try:
            info_dict[key_val[0]] = eval(key_val[1])
        except:
            info_dict[key_val[0]] = key_val[1]

    return info_dict

info_dict = create_info_dict(df1['INFO'][0])
info_dict


######## Main Function ############
def variant_annotator_main():
    """Reads in a VCF file and writes out an annotated variant file."""

    return -1

