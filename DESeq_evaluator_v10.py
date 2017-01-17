import csv
import argparse
import math
from collections import defaultdict

def makeheader(inlist):
    header =  ['gene',' | ']
    seperator=['********']
    for element in inlist:
        header.append(element.name)
        header.append('|')

        for r in range(len(element.name)):
            seperator.append('*')

    sep = ''.join(seperator)
    return header, sep

def contains(list_, filter_):
    """
Find object in List of objects
    """
    for x in list_:
        if filter_(x):
            return True
    return False

def fileparser(string):
    genelist = []
    with open((string)) as in_raw:
        infile = csv.reader(in_raw)
        infile.next()
        for row in infile:
            new_gene = gene(row[0])
            new_gene.logfold.append(float(row[2]))
            new_gene.basemean.append(float(row[1]))
            genelist.append(new_gene)
    return genelist


def getlabels(infile_list):
    labels = []
    for inf in infile_list:
        labels.append(inf.name)
    return labels
    
def getratios(results_condition,conditions):
    """
Returns sets of DE genes that are active in certain conditions
    """
    setlist = {}
    for r in range(conditions):
        setlist[r] = []
    
    for gene in results_condition.genelist:
        conditions = len(gene.logfold)

        count = 0
        for set_ in setlist:
            if gene.logfold[count] > 1.5 or gene.logfold[count] < -1.5 :
                setlist[count].append(gene.name)
            count +=1
    return setlist
        
        
        
    

class gene():
    """
Genes hava an unique ID, a list of log2foldchanges, and a list of basemeans
    """

    def __init__(self, name):
        self.name = name
        self.logfold = []
        self.basemean = []


class condition():
    """
conditions are DESeq2 outputs,  storing them here to get a stable comparison
conditions contain a name and a list of genes with logfold and expression
output can be a condition
    """

    
    def __init__(self, name, genelist ):
        self.name = name
        self.genelist = genelist
        self.header = []


parser = argparse.ArgumentParser(description='Parses DESeq 2 output')
parser.add_argument('-i','--infile', dest='infiles', type=str, nargs='+',
                    help='One or more DESeq input files')
parser.add_argument('-o','--output', dest='outfile', type=str,
                    default='none',
                    help='[OPTIONAL] writers output to file,  expects name for output file,  default= do not write to file' )

parser.add_argument('-bm','--basemean', dest='basemean', type=int,
                    default=100,
                    help='[OPTIONAL] Cutoff for basemean,  default 100' )

parser.add_argument('-l','--logfoldchange', dest='log2fold', type=float,
                    default= 1.5,
                    help='[OPTIONAL] Cutoff for log2foldchange,  default +-1.5' )
parser.add_argument('-d','--direction', dest='direction', type=str,
                    default='both',
                    help="[OPTIONAL] Looking for 'up' / 'down' regulated genes,  or default 'both'" )
parser.add_argument('-v','--venndiagram', dest='venn', type=bool,
                    default=False,
                    help="[OPTIONAL] Want to plot a simple venn diagram ? " )


##parser.add_argument('-c','--conversion', dest='convers', type=file,
##                    const=sum, default='none,
##                    help='[OPTIONAL] conversion table to get Ortholog gene names ')
##
##

args = parser.parse_args()
print(args)


infile_list = []

for infile in args.infiles:
    gene_list = fileparser(infile)
    infile_list.append(condition(infile, gene_list ))


# initiate some default behaviour
if len(infile_list)== 2:
    getdelta = True
else:
    getdelta =  False

# parse arguments
writeout = False
if args.outfile != 'none':
    writeout = True
    


    

output = condition('output',[])

logfiltering = args.direction

# go through the individual files get genes of interest and put them into a new outputfile

for element in infile_list:
    #print len(element.genelist), element.name
    
    for individual_gene in element.genelist:
        try:
            if logfiltering == 'both':
                express = math.sqrt((float(individual_gene.logfold[0])*float(individual_gene.logfold[0])))
            if logfiltering == 'up':
                express = (float(individual_gene.logfold[0]))
            if logfiltering == 'down':
                express = (float(individual_gene.logfold[0]) * -1)
            else:
                #print ' no valid direction provided,  printing both up and down regulated genes '
                express = math.sqrt((float(individual_gene.logfold[0])*float(individual_gene.logfold[0])))
        except:
            express = 0

        if individual_gene.basemean[0] > args.basemean and express > args.log2fold:
            # add empty genes to the output gene list, populate values later
            if not contains( output.genelist, lambda gene: gene.name == individual_gene.name):
                output.genelist.append(gene(individual_gene.name))



# now populate the genes of interest with values from the files
for element in output.genelist:
    for infile in infile_list:
        element.logfold.append(filter(lambda gene: gene.name == element.name, infile.genelist)[0].logfold[0])
        

header, seperator = makeheader(infile_list)
print seperator
print ' '.join(header)
print seperator

        
for element in output.genelist:
    outlist = [element.name]
    for logfol in element.logfold:
        outlist.append(str(logfol))
    print  ' | '.join(outlist)
                       

if writeout:
    with open('%s'%(args.outfile)) as out_raw:
        outfile = csv.writer(out_raw)

        head = []
        for element in header:
            if not ' | ' in element:
                head.append(element)
        outfile.writerow(head)

        for element in output.genelist:
            outlist = [element.name]
            for logfol in element.logfold:
                outlist.append(str(logfol))
            outfile.writerow(outlist)
              

#  plot the venn diagram if desired:
if args.venn :
    print 'Creating Exploratory Venn diagram from the data'
    from matplotlib import pyplot as plt
    import numpy as np

    plt.figure(figsize=(10,10))

    plt.title("Exploratory venn diagram")
    labels = getlabels(infile_list)
    ratios = getratios(output, len(infile_list))
##    for key, value in ratios.items():
##        print 'ratio'  , key, value
    try:
        if len(infile_list) == 2:
            from matplotlib_venn import venn2,  venn2_circles
            set1 = set(ratios[0])
            set2 = set(ratios[1])
            venn2([set1, set2], (infile_list[0].name,infile_list[1].name))
            c = venn2_circles(subsets=[set1, set2], linestyle='dashed')
            plt.show()
        if len(infile_list) == 3:
            from matplotlib_venn import venn3,  venn3_circles
            set1 = set(ratios[0])
            set2 = set(ratios[1])
            set3 = set(ratios[2])
            venn3([set1, set2, set3], (infile_list[0].name,infile_list[1].name,infile_list[2].name))
            c = venn3_circles(subsets=[set1, set2, set3], linestyle='dashed')
            plt.show()
    except:
        print 'Matplotlib venn not found,  please install using the following command :'
        print 'easy_install matplotlib-venn'
        exit

    

    
    
















