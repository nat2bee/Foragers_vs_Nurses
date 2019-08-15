#!/usr/local/bin/python


"""
Take an annotation output from Annocript and print it in the format used as input in TopGo custom annotation ".map".
To be read in R using: geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))

Usage = Annocript2TopGo.py -c <ontology_category> -a <annotation> -o <output>

Where: 
ontology_category = The ontology category you desire as annotation (BP: Biological Process; MF: Molecular function; CC: Cellular Component or all) 
                    If more then one category is desired you can use two symbols separated by comma ex. "BP,MF". For all terms use "all"
annotation = table result from Annocript containing the information from annotation (...filt_ann_out.txt)
output = the name of the output to save the table with the information of the transcripts of interest

Options: -h for usage help
"""

import sys, getopt

# Useful variables

n = 0
BP= False
MF= False
CC= False
all= False



# Check for the arguments, open the inputs and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"hc:a:o:",["annotation","output"])
except getopt.GetoptError:
    print '\n', '####     Invalid use     ####', '\n'
    print 'Usage = Annocript2TopGo.py -c <ontology_category> -a <annotation> -o <output>'
    print 'For help use Annocript2TopGo.py -h'
    sys.exit(99)
    
for opt, arg in opts:
    if opt == '-h':
        print '\n', 'Take a annotation output from annocript and print it in the format used as input in TopGo for custom annotation ".map".'
        print 'To be read in R using: geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))', '\n'
        print 'Usage = Annocript2TopGo.py -c <ontology_category> -a <annotation> -o <output>'
        print 'Where:'
        print 'ontology_category = The ontology category you desire as annotation (BP: Biological Process; MF: Molecular function; CC: Cellular Component or all)'
        print 'If more then one category is desired you can use two symbols separated by comma ex. "BP,MF". For all terms use "all"'
        print 'annotation = table result from Annocript containing the information from annotation (...filt_ann_out.txt)'
        print 'output = the name of the output to save the table with the information of the transcripts of interest'
        sys.exit()
    if opt in ("-c"):
        category = arg
        if category in "all":
            all = True
        elif "," in category:
            category = category.split(",")
        elif "BP" in category:
            BP = True
        elif "MF" in category:
            MF = True
        elif "CC" in category:
            CC = True
    if opt in ("-a", "--annotation"):
        annotation = open(arg)
    if opt in ("-o", "--output"):
        outname = arg + ".map"
        output = open(outname,"w")

  



## Open the file and save the gene id and the ontology category of interest in the output file

for line in annotation: 
    if n == 0: # skip the header
        n = 1
        continue
    line = line.replace("]---[", ",")  # Change the separator
    elements = line.split("\t")
    
    transcript_id = str(elements[0]) # Interesting data
    GO_BP = str(elements[27])
    GO_BP = GO_BP.replace("-", " ")  # Change blank annotation symbol 
    GO_MF = str(elements[29])
    GO_MF = GO_MF.replace("-", " ")  
    GO_CC = str(elements[31])
    GO_CC = GO_CC.replace("-", " ")  

    
    # save only what the user want
    if all is True:
        all_GO = GO_BP + "," + GO_MF + "," + GO_CC  
    elif BP is True and MF is True and CC is False:
        all_GO = GO_BP + "," + GO_MF 
    elif BP is True and CC is True and MF is False:
        all_GO = GO_BP + "," + GO_CC
    elif MF is True and CC is True and BP is False:
        all_GO = GO_MF + "," + GO_CC
    elif BP is True and CC is False and MF is False:
        all_GO = GO_BP 
    elif CC is True and BP is False and MF is False:
        all_GO = GO_CC
    elif MF is True and CC is False and BP is False:
        all_GO = GO_MF
    output.write(transcript_id + "\t" + all_GO + "\n" )
    
    
output.close()
