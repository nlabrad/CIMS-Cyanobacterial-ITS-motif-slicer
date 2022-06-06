import re, sys, argparse
from Bio import Entrez, __version__, SeqIO #import the required libraries

#TODO: add logging.
#TODO: handle not having an argument after the option -f or -g
#TODO: handle not finding the file passed.

# Gets a sequence from genbank
def fetch_genbank(accession): #fetch sequence and taxonomy with accession#
    Entrez.email = 'nlab@fastmail.com' # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
        seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        fasta = (">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))
        print(fasta)
    return fasta #pass this to parse_fasta to get sequence


# Gets sequences from a fasta file. Calls findMotif for each.
def parse_fasta(fasta):
    for seq_record in SeqIO.parse(fasta, "fasta"): #for each entry do the below.
        print("Organism name: ", seq_record.id)
        its_sequence = extractITS(str(seq_record.seq))
        if(extractITS == "None"):
            print("No ITS region found in this sequence")
        else:
            # print("\nITS Region sequence: ", its_sequence)
            print("\nITS Region length: ", len(its_sequence))
            motifs = findMotifs(its_sequence)
            print("Motif  |   [ Start, End, Sequence, Lenght] \n")
            for key in motifs:
                if(motifs[key] == None):
                    print(key, ": Not found\n")
                else:
                    print(key, ': ', motifs[key],'\n')
            print('------------------------------------------------------------------------------------\n')

#Extraction of ITS region from a 16S-23S region.
# Find CCTCCTT, get rid of what's before that.
def extractITS(geneRegion):
    if(re.search("CCTCCTT", geneRegion) == None):
        print("ITS Region not found in this sequence")
        return None
    else:    
        itsStartPosition = re.search("CCTCCTT", geneRegion).start() + 7
        itsseq = geneRegion[itsStartPosition:]
        return itsseq

#Motifs reference:
# Motifs["key"] = {'ITS': its region
#                  'motifName': [start_position, end_position, sequence_string, length]}

def findMotifs(its_region): #Find the motifs
    motifs = {} #dictionary. motifs[motif-name]=[start-position,end-position, sequence, length]
    
    motifs["ITS"] = its_region #store the whole ITS region to be used as a reference.
    index_shift = 0
    
    d1d1Search = re.search(r"GACCT(.*?)AGGTC", its_region) #find text between basal clamps, starting with GACCT/C to the first AGGTC (*? is lazy search)
    if (d1d1Search == None): 
        d1d1Search = re.search(r"GACCA(.*?)[AT]GGTC",its_region) 
    if (d1d1Search == None): 
        d1d1Search = re.search(r"GACCG(.*?)CGGTC",its_region)
    if (d1d1Search == None):
            d1d1Search = None
    else: 
        motifs["leader"] = [0,d1d1Search.start()] #start-end positions relative to the ITS region
        motifs["leader"].append(motifs["ITS"][motifs["leader"][0]:motifs["leader"][1]]) #The sequence string of the motif
        motifs["leader"].append(len(motifs["leader"][2])) #The length of the motif
        
        motifs["d1d1"] = [d1d1Search.start(),d1d1Search.end()] #[0] = start, [1] = end
        motifs["d1d1"].append(motifs["ITS"][motifs["d1d1"][0]:motifs["d1d1"][1]]) #[2] = sequence string
        motifs["d1d1"].append(len(motifs["d1d1"][2])) #[3] = lenght
        
        index_shift = d1d1Search.end() #count of processed characters. Used to keep track of absolute position in the index.
        its_region = its_region[d1d1Search.end():] #trim processed its region to avoid overlapping searches.
    
    
    # spacer-D2-spacer D3-spacer regionand  tRNA-Ile(1)
    tRNA1Search = re.search(r"GGGC[TC]ATTA(.*?)GGCCCA",its_region) #find text between basal clamps
    if (tRNA1Search == None):
        motifs["tRNA1"] = None
    else:  
        motifs["sp_d2d3_sp"] = [d1d1Search.end(), tRNA1Search.start()+index_shift]
        motifs["sp_d2d3_sp"].append(motifs["ITS"][motifs["sp_d2d3_sp"][0]:motifs["sp_d2d3_sp"][1]])
        motifs["sp_d2d3_sp"].append(len(motifs["sp_d2d3_sp"][2]))
        
        motifs["tRNA1"] = [tRNA1Search.start()+index_shift, tRNA1Search.end()+index_shift]
        motifs["tRNA1"].append(motifs["ITS"][motifs["tRNA1"][0]:motifs["tRNA1"][1]])
        motifs["tRNA1"].append(len(motifs["tRNA1"][2]))
                        
        index_shift = index_shift + tRNA1Search.end() #add to the index_shift what will be trimmed in the next line
        its_region = its_region[tRNA1Search.end():] #trim processed its region
    
    #TODO: V2 region is between these
    
    #find tRNA-Ala(2)
    tRNA2Search = re.search(r"GGGG(.*?)[TC]CTCCA",its_region)  #find text between basal clamps
    if (tRNA2Search == None):
        motifs["tRNA2"] = None
    else: 
        motifs["tRNA2"] = [tRNA2Search.start()+index_shift, tRNA2Search.end()+index_shift]
        motifs["tRNA2"].append(motifs["ITS"][motifs["tRNA2"][0]:motifs["tRNA2"][1]])
        motifs["tRNA2"].append(len(motifs["tRNA2"][2])) 
        
        index_shift = index_shift + tRNA2Search.end()
        its_region = its_region[tRNA2Search.end():] #trim processed its region
        
    #find BoxB
    BoxBSearch = re.search(r"[TC]AGCA[AC](.*?)TGCT[AG]", its_region) #find text between basal clamps
    if (BoxBSearch == None):
        motifs["BoxB"] = None
    else:        
        motifs["BoxB"] = [BoxBSearch.start()+index_shift, BoxBSearch.end()+index_shift]
        motifs["BoxB"].append(motifs["ITS"][motifs["BoxB"][0]:motifs["BoxB"][1]])
        motifs["BoxB"].append(len(motifs["BoxB"][2]))
        index_shift = index_shift + BoxBSearch.end()
        its_region = its_region[BoxBSearch.end():] #trim remaining its region
        
    #find Box-A
    BoxASearch = re.search(r"G[CA]A(.*?)G[AG]AAACT",its_region) #TODO: make it search backwards from ACT and find the first pattern matched.
    if (BoxASearch == None):
        motifs["BoxA"] = None
    else:  
        motifs["BoxA"] = [BoxASearch.start()+index_shift, BoxASearch.end()- 3 +index_shift] #-3 to exclude ACT
        motifs["BoxA"].append(motifs["ITS"][motifs["BoxA"][0]:motifs["BoxA"][1]])
        motifs["BoxA"].append(len(motifs["BoxA"][2]))
        index_shift = index_shift + BoxASearch.end() - 3
        its_region = its_region[BoxASearch.end() - 3:] #trim remaining its region   
    
    
    #find D4  
    D4Search = re.search(r"ACT(.*?)TA[TGAC]",its_region) #Should be right after the GAA above
    if (D4Search == None):
        motifs["D4"] = None
    else:  
        motifs["D4"] = [D4Search.start()+index_shift, D4Search.end()+index_shift]
        motifs["D4"].append(motifs["ITS"][motifs["D4"][0]:motifs["D4"][1]])
        motifs["D4"].append(len(motifs["D4"][2]))
        index_shift = index_shift + D4Search.end()
        its_region = its_region[D4Search.end():] #trim remaining its region
        
    #Find v3
    V3Search = re.search(r"GT[CA]G(.*?)CA[CG]A[GC]",its_region) 
    if (V3Search == None):
        motifs["v3"] = None
    else:  
        motifs["v3"] = [V3Search.start()+index_shift, V3Search.end()+index_shift]
        motifs["v3"].append(motifs["ITS"][motifs["v3"][0]:motifs["v3"][1]])
        motifs["v3"].append(len(motifs["v3"][2]))
        index_shift = index_shift + V3Search.end()
        its_region = its_region[V3Search.end():]
        
    return motifs
    # end of extractMotif
   
   
# If -g then run fetch then parse_fasta, 
# If -f then just run parse_fasta

# Create the parser    
parser = argparse.ArgumentParser(description = "Find motifs within Cyanobacteria ITS region sequences")

#Create a group of mutually exclusive argument options
group = parser.add_mutually_exclusive_group()

#Add arguments to group
group.add_argument('-f', '--fasta', action="store", help = "Find motifs in a given fasta file.")
group.add_argument('-g', '--genbank', help = "Fetch a sequence from a given genbank accession number.")

#Add standalone argument
parser.add_argument('-H', '--html', help = "Create an HTML file with the motifs highlighted over the whole sequence.")

#If there is no argument then print the help
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

#Parse the arguments, store them in args.
args = parser.parse_args()
    
if(args.fasta):
    #TODO: make sure the file exists.
    try:
        fasta = open(args.fasta, "r")
        print("Processing fasta file...")
        parse_fasta(fasta)
        exit()
    except IOError:
        print ("File not found")
        exit() 
    
if(args.genbank):
    print("Fetching genbank data...")
    fasta = fetch_genbank(args.genbank)
    parse_fasta(fasta)
