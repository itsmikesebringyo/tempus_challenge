# Mike Sebring
# Tempus Coding Challenge
# October 19, 2020

import requests
import pandas as pd


filepath = 'Challenge_data_(1).vcf'


## Helper function.. read vcf file to pandas dataframe.
def vcf_to_dataframe(filepath):
    """Takes a vcf filepath as input and returns a variant dataframe with a Variant ID column."""
    
    # Open the file and find the starting point of the variant data.
    # Set the cursor position and read in as a dataframe.
    with open(filepath) as f:
        start = f.read().find('#CHROM')
        f.seek(start+1)
        df1 = pd.read_csv(f, sep='\t')

    # Next, create the variant ID column which we'll use to annotate.
    df1['variant_id'] = df1.agg(lambda x: f"{x['CHROM']}-{x['POS']}-{x['REF']}-{x['ALT']}", axis=1)

    return df1

# testing dataframe
variant_dataframe = vcf_to_dataframe(filepath)

variant_dataframe.columns
variant_dataframe
variant_dataframe['CHROM'][2]

variant_dataframe['variant_id'][1000]

## Helper function to get allele frequency from ExAC API.

def get_allele_freq(variant_id, vcf_AF):
    """Takes variant ID and AF value as input and returns allele frequency."""
    
    # Hit the ExAC API and convert the results into JSON.
    variant_info = requests.get(f'http://exac.hms.harvard.edu/rest/variant/variant/{variant_id}').json()
    
    # Try getting the allele frequency from the API result if possible.
    # Otherwise, we need to use the allele frequency provided in source file.
    try:
        allele_frequency = variant_info['allele_freq']
    except:
        allele_frequency = vcf_AF

    return allele_frequency

# testing ExAC API
test_allele_freq = get_allele_freq(variant_dataframe['variant_id'][1000], 0.5)
print(test_allele_freq)

## Helper function to create a dictionary of values from the INFO column.

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

# testing INFO dictionary
info_dict = create_info_dict(variant_dataframe['INFO'][0])
info_dict

## Helper function to create output file.

def dataframe_to_outfile(df1, outfile_name, delim='\t'):
    """Takes dataframe, filename, and delimiter as input and creates file."""

    df1.to_csv(outfile_name, sep=delim)

######## Main Function ############
def variant_annotator_main(vcf_path):
    """Reads in a VCF file and writes out an annotated variant file."""

    # First, we need to read the vcf file and create the dataframe.
    variant_df = vcf_to_dataframe(vcf_path)

    # Now we need to iterate through each variant in the file, grab the 
    # necessary annotations, and append them to a new dataframe.
    
    # Create new dataframe.
    col_list = ['variant_id', 'variant_type', 'variant_effect', 'coverage_depth', 'variant_reads', 
                'variant_percent', 'allele_frequency', 'additional_info']

    output_df = pd.DataFrame(columns=col_list)

    for i in range(len(variant_df)):

        # Create a dictionary out of the INFO column so the values
        # can be used to annotate the current variant.
        info_dict = create_info_dict(variant_df['INFO'][i]) 

        # Variant type is found with the TYPE attribute
        variant_type = info_dict['TYPE']
        # Variant effect TBD
        variant_effect = -1
        # Depth of sequence coverage is found with the DP attribute
        coverage_depth = info_dict['DP']
        # Number of reads supporting variant is found with addition of SAF and SAR attributes
        variant_reads = info_dict['SAF'] + info_dict['SAR']
        # The percentage of reads supporting variant vs those supporting reference reads
        # can be found using the following calculation of attributes: (SAF + SAR) / (SAF + SAR + SRF + SRR)
        variant_percent = (variant_reads) / (variant_reads + info_dict['SRF'] + info_dict['SRR'])
        variant_percent = float(f"{variant_percent:.2f}")
        # Get the allele frequency using our helper function. Use the estimated 
        # alle frequency attribute (AF) as backup
        allele_freq = get_allele_freq(variant_df['variant_id'][i], info_dict['AF'])

        # Create a new row to be appended to the output dataframe.
        new_row = {'variant_id': variant_df['variant_id'][i],
                   'variant_type': variant_type,
                   'variant_effect': variant_effect,
                   'coverage_depth': coverage_depth,
                   'variant_reads': variant_reads,
                   'variant_percent': variant_percent,
                   'allele_frequency': allele_freq,
                   'additional_info': -1}

        # Append row to our output dataframe
        output_df = output_df.append(new_row, ignore_index=True)
        
    dataframe_to_outfile(output_df, 'output_file.tsv')

    # End of main.



variant_annotator_main(filepath)