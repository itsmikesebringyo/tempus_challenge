# Mike Sebring
# Tempus Coding Challenge
# October 22, 2020

# Import necessary modules.
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

    return df1


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

        info_dict[key_val[0]] = key_val[1]

    return info_dict


## Helper function to create output file.

def dataframe_to_outfile(df1, outfile_name, delim='\t'):
    """Takes dataframe, filename, and delimiter as input and creates file."""

    df1.to_csv(outfile_name, sep=delim)


## Helper function to split attributes in the instance where multiple types are given.

def multiple_type_split(info_attribute):
    """Takes variant attribute as input and returns first value."""

    # Split the inputted attribute and return the first value.
    attirbute_split = info_attribute.split(',')

    return attirbute_split[0]


## Helper function to get the most deleterious effect based on variant type.

def get_variant_effect(variant_type):
    """Takes variant type as input and returns the most deleterious effect."""

    # Set up a mapping dictionary of variant types and their most deleterious effect.
    mapping = {
        'snp': 'nonsense',
        'mnp': 'nonsense',
        'ins': 'frameshift',
        'del': 'frameshift',
        'complex': 'frameshift'
    }

    # Get the variant effect and return the value.
    variant_effect = mapping[variant_type]

    return variant_effect


##################################### Main Function ####################################

def variant_annotator_main(vcf_path):
    """Reads in a VCF file and writes out an annotated variant file."""

    # First, we need to read the vcf file and create the dataframe.
    variant_df = vcf_to_dataframe(vcf_path)

    # Now we need to iterate through each variant in the file, grab the 
    # necessary annotations, and append them to a new dataframe.
    
    # Create output dataframe.
    col_list = ['variant_id', 'variant_type', 'variant_effect', 'coverage_depth', 'variant_reads', 
                'variant_percent', 'allele_frequency', 'additional_info']

    output_df = pd.DataFrame(columns=col_list)

    for i in range(250): ########## substitute in the length of the dataframe when finalized

        # Create a dictionary out of the INFO column so the values
        # can be used to annotate the current variant.
        info_dict = create_info_dict(variant_df['INFO'][i]) 

        # Next, we need to run a check to see if there are multiple types.
        multiple_vals = False
        if len(info_dict['TYPE'].split(',')) > 1:
            multiple_vals = True

        # Now create variant_id which will be used to annotate and call ExAC API.
        # Use helper function to return first ALT allele if there are multiple.
        if multiple_vals:
            variant_id = "-".join([str(variant_df['CHROM'][i]), str(variant_df['POS'][i]), variant_df['REF'][i], multiple_type_split(variant_df['ALT'][i])])
        else:
            variant_id = "-".join([str(variant_df['CHROM'][i]), str(variant_df['POS'][i]), variant_df['REF'][i], variant_df['ALT'][i]])

        # Variant type is found with the TYPE attribute. 
        # Use helper function to return first value if there are multiple.
        if multiple_vals:
            variant_type = multiple_type_split(info_dict['TYPE'])
        else:
            variant_type = info_dict['TYPE']
        
        # Get the variant effect based on the type.
        variant_effect = get_variant_effect(variant_type)
        
        # Depth of sequence coverage is found with the DP attribute
        coverage_depth = eval(info_dict['DP'])
        
        # Number of reads supporting variant is found with addition of SAF and SAR attributes.
        # Use helper function to return the first value if there are multiple.
        if multiple_vals:
            variant_reads = eval(multiple_type_split(info_dict['SAF'])) + eval(multiple_type_split(info_dict['SAR']))
        else:
            variant_reads = eval(info_dict['SAF']) + eval(info_dict['SAR'])
        
        # The percentage of reads supporting variant vs those supporting reference reads
        # can be found using the following calculation of attributes: (SAF + SAR) / (SAF + SAR + SRF + SRR)
        variant_percent = (variant_reads) / (variant_reads + eval(info_dict['SRF']) + eval(info_dict['SRR']))
        variant_percent = float(f"{variant_percent:.2f}")
        
        # Get the allele frequency using our helper function. Use the estimated 
        # allele frequency attribute (AF) as backup
        if multiple_vals:
            backup_AF = eval(multiple_type_split(info_dict['AF']))
        else:
            backup_AF = eval(info_dict['AF'])
        
        allele_freq = get_allele_freq(variant_id, backup_AF)

        # For additional information, let end consumer know the variant had multiple types.
        if multiple_vals:
            additional_info = 'multiple types'
        else:
            additional_info = ''

        # Create a new row to be appended to the output dataframe.
        new_row = {'variant_id': variant_id,
                   'variant_type': variant_type,
                   'variant_effect': variant_effect,
                   'coverage_depth': coverage_depth,
                   'variant_reads': variant_reads,
                   'variant_percent': variant_percent,
                   'allele_frequency': allele_freq,
                   'additional_info': additional_info}

        # Append row to our output dataframe
        output_df = output_df.append(new_row, ignore_index=True)

        # End of for-loop.
        
    dataframe_to_outfile(output_df, 'output_file.tsv')

    # End of main.


variant_annotator_main(filepath)