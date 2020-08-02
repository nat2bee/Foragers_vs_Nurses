#!/usr/local/bin/python


"""
## 17 July 2020 ##

Filter orthogroups ( Orthogroups.tsv ) containing a set of genes of interest (gene_list).

Developed (in my case) to use with the result table from OrthoFinder "Orthogroups.tsv" 
to get the orthogroups containing the differentially expressed genes.

Usage = filter_orthogroups.py <options> -s <species_name> -i <gene_list> -t <Orthogroups.tsv> 

Where: 
species_name = name of the species from which your gene_list is from, as it appears in the Orthogroups.tsv column
gene_list = list with all the genes name (one per line) you are interested
Orthogroups.tsv = table output from Orthofinder (Should be in the folder the script will be run)

Options: 
-h for usage help
-n to indicate the number of proteins in the protein set (not the same as the number of genes)


Output: filtered table will be saved under the name <species_name_filtered_Orthogroups.tsv>
"""

import sys, getopt

sps_name = ""
n_used = False

# Check for the arguments, open the inputs and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"hN:s:i:t:",["list=","table="])
except getopt.GetoptError:
    print '\n', '####     Invalid use     ####', '\n'
    print 'Usage = filter_orthogroups.py <options> -s <species_name> -i <gene_list> -t <Orthogroups.tsv>'
    print 'For help use filter_orthogroups.py -h'
    sys.exit(99)
    
for opt, arg in opts:
    if len(arg) < 3 and opt == '-h':
        print '\n', 'Filter orthogroups ( Orthogroups.tsv ) containing a set of genes of interest (gene_list).', '\n'
        print 'Usage = filter_orthogroups.py <options> -s <species_name> -i <gene_list> -t <Orthogroups.tsv>'
        print 'Where: species_name = name of the species from which your gene_list is from, as it appears in the Orthogroups.tsv column'
        print 'gene_list = list with all the genes name (one per line) you are interested'
        print 'Orthogroups.tsv = table output from Orthofinder\n'
        print 'Options: -h for usage help'
        sys.exit()
    elif len(arg) >= 3:
    	if opt in ("-N"):
    		n_used = True
    		n_prot = int(arg)
    	if opt in ("-s"):
    		sps_name = str(arg)
        if opt in ("-i", "--list"):
            input_list = open(arg)
        if opt in ("-t", "--table"):
            ortho_table = open(arg)
            out_name = sps_name + "_filtered_" + arg
    elif len(arg) < 3:
        print '\n', '###    Arguments are missing   ###', '\n', '\n' 'Use -h option for help\n'
        sys.exit(1)
    else:
        assert False, "unhandled option"

## Open output file
output = open(out_name,"w")

## other variables
header = True
gene_name = ""    
gene_list = []
ortho_prots = []
ortho_genes = []
n = 0
    

# Make a list with input gene names
for name in input_list:
    gene_name = name.split("\n")
    gene_name = str(gene_name[0]) + "." # to make it more specific
    gene_list.append(gene_name)


# Check in each row of the table if this name appears in the species column
for line in ortho_table:
	if header is True: 
		output.write(line) # save the header
		split_line = line.split("\r")
		split_line = split_line[0]
		split_line = split_line.split("\t")
		column = split_line.index(sps_name) # get the location of the sps in the file
		header = False
		continue
	
	else: 
		split_line = line.split("\r")
		split_line = split_line[0]
		split_line = split_line.split("\t")
		genes_orthos = split_line[column]
		genes_orthos = genes_orthos.split(", ")
		new_line = True # mark it is a new orthogroup 
		
		# To each protein in the orthogroup
    	for ortho in genes_orthos:
    		# get the gene name
    		new_ortho = ortho.split("p")
    		new_ortho = new_ortho[0]
    		
    		# Check if the gene is in the list of interest
    		if new_ortho in gene_list:
    			
    			# if it is save the orthogroup, only once
    			if new_line is True:
    				output.write(line)
    				new_line = False # mark the orthogroup is saved
    			
    			# If this protein appears for the first time, register it
    			if ortho not in ortho_prots:
    				ortho_prots.append(ortho)
    				
    			# If the gene appears for the first time, register it
    			if new_ortho not in ortho_genes:
    				ortho_genes.append(new_ortho)


# remove possible repetitions of genes
ortho_genes_uniq = set(ortho_genes)	

print "\n>>>\nNumber of genes in the list: ",   len(gene_list) 
print "Number of proteins in Orthogroups: ", len(ortho_prots)
print "Number of genes in Orthogroups: ", len(ortho_genes_uniq)
print "% of genes in orthogroups: ", float(len(ortho_genes_uniq)*100)/len(gene_list) , "%\n>>>\n"

if n_used is True:
	print "Number of proteins in the list: ",   n_prot
	print "% of proteins in orthogroups: ", float(len(ortho_prots)*100)/n_prot , "%\n>>>\n"

    		

output.close()
