# Tempus Coding Challenge

The intent of this script is to ingest a VCF file and create an output file that annotates each variant from the input file. Each variant is annotated with the following information:

1. The type of variation
2. The variation effect (if multiple types, the most deleterious effect was used)
3. The depth of sequence coverage at the site of variation
4. The number of reads supporting the variant
5. The percent of reads supporting the variant versus those supporting reference reads
6. The allele frequency of the variant from the ExAC API
7. An additional annotation stating whether there were multiple types

## Assumptions

Some assumptions were made during the creation of this script.

1. The first assumption is that the variable type is found in the INFO column of the VCF file after the TYPE attribute.
2. The second assumption is that the variable effect depends on the variable type. Variants with type "snp" and "mnp" mapped to the most deleterious effect of a "nonsense" mutation. Variants with type "ins", "del", and "complex" mapped to the most deleterious effect of a "frameshift" mutation.
3. The third assumption is that the depth of sequence coverage is found in the INFO column of the VCF file after the DP attribute.
4. The fourth assumption is that the number of reads supporting the variant is found in the INFO column of the VCF file and is the addition of the SAF and SAR attributes.
5. The fifth assumption is that the percentage of reads supporting the variant versus those supporting reference reads is found in the INFO column of the VCF file and is the calculation of the following attributes: (SAF + SAR) / (SAF + SAR + SRF + SRR)
6. The sixth assumption is that the first TYPE of each variant is the most deleterious.

## Details

The script is composed of a main function along with several helper functions.

The main function starts by calling a helper function that reads in the VCF file into a pandas dataframe. Next, an empty dataframe is created (with the necessary annotations as columns) as a placeholder so that we can add rows to it.

Next, we start iterating through each row of the VCF dataframe and getting the annotations as we go. As we iterate through each variant, a dictionary is created from the INFO column since we are using several of the attributes from that column. Once we've created our dictionary, we check to see if there are multiple types for the current variant in the iteration. We use this knowledge along with the sixth assumption from above; if there are multiple types for the variant, we pull the first value of each attribute where multiple values are provided.

The next step was to create the "variant_id". This served two purposes:

It serves as the first column in the output file
We use it to call the ExAC API
From here, several more helper functions are called to get the variables needed for annotation. Once these variables are gathered, they are put into a dictionary and this dictionary is added as a new row to the output dataframe.

After the input VCF dataframe has been completely iterated through and the annotations have been written to the output dataframe, the output dataframe is then written to an output file for the end consumer.

## References

The following references helped created the assumptions from above and provided insight into the formatting of VCF files:

* http://www2.csudh.edu/nsturm/CHEMXL153/DNAMutationRepair.htm
* https://genome.sph.umich.edu/wiki/Variant_classification
* https://samtools.github.io/hts-specs/VCFv4.1.pdf
* http://www.htslib.org/doc/vcf.html
* http://exac.hms.harvard.edu/

