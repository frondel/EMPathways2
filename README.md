# EMPathways2
EMPathways2 - A tool for pathway reconstruction from metagenomic data

Requirements
Python 3.x

# Usage
python EMPathways2.py [-h] -ge GENE_ENZYME_FILE -epd ENZYME_PATHWAY_DICTIONARY -ge GENE_EXPRESSION_FILE [-eo ENZYME_OUTPUT] [-po PATHWAY_OUTPUT] [--theta_e THETA_E] [--theta_p THETA_P]

Arguments
-h, --help : Show the help message and exit.

-ge, --gene_enzyme_file : Required. The gene-enzyme relation file (gene-enzyme-weight file with 1's).

-epd, --enzyme_pathway_dictionary : Required. The enzyme pathway dictionary (EPD) file.

-ge, --gene_expression_file : Required. The gene expression file (gene_fpkm IsoEM output file).

-eo, --enzyme_output : Optional. The filename to write enzyme data to. If not provided, the output will be printed to the console.

-po, --pathway_output : Optional. The filename to write pathway data to. If not provided, the output will be printed to the console.

--theta_e : Optional. Theta value for enzyme convergence. Default is 0.000001.

--theta_p : Optional. Theta value for pathway convergence. Default is 0.000001.


# Example
python EMPathways2.py -ge gene_enzyme_file.txt -epd enzyme_pathway_dictionary.txt -ge gene_expression_file.txt -eo enzyme_output.txt -po pathway_output.txt --theta_e 0.00001 --theta_p 0.00001
