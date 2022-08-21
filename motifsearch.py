import os
import re, sys, argparse, colorama, json
from Bio import Entrez, SeqIO #import the required libraries
from colorama import Fore, Back, Style
colorama.init(autoreset=True)

#TODO: add logging.
#TODO: option to create CSV
#TODO: HTML with highlighted motif text
#TODO: option to select just one motif from given seq instead of all the motifs.
#TODO: option to exclude ones without tRNAs, or alternatively, just include those with both tRNAs

def generate_json(dict):
    jsonString = json.dumps(motifs, indent=3)
    file = open("motifs.json", "w")
    file.write(jsonString)
    file.close()
    
def print_motifs(motifs):
        if(motifs == None):
            print(Back.RED + Fore.WHITE + "ITS Region not found in this sequence!")
        else:
            for key in motifs:
                if(motifs[key] == None):
                    if (key == "tRNA1" or key == "tRNA2"):
                        print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not present in this operon.")
                        print("\n")  
                    else:  
                        print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not found in this sequence.")
                        print("\n")
                else:
                    #if (key == 'd1d1' or key == 'BoxB'):
                    if len(motifs[key]) > 1:
                        print(Fore.CYAN + Style.BRIGHT + key + " " +Fore.RED + Style.BRIGHT + str(len(motifs[key])) + " possible sequences found: \n")     
                        for index, item in enumerate(motifs[key]):
                            print(Fore.MAGENTA + Style.BRIGHT +  "Sequence " + str(index+1) + ": " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                                Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")
                    else:
                        print(Fore.CYAN + Style.BRIGHT + key + ":")
                        for item in motifs[key]:
                            print(Fore.MAGENTA + Style.BRIGHT +  "Sequence " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                                Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")
                    # else:
                    #     print(Fore.CYAN + Style.BRIGHT + key + "\n" + 
                    #         Fore.MAGENTA +  "Sequence: " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(motifs[key]) + 
                    #         Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(motifs[key])) + "\n")    

# Gets the specified sequence from genbank to process.
def parse_genbank(accession, email): #fetch sequence and taxonomy with accession#
    Entrez.email = email # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
        seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        motifs = slice_motifs(str(seq_record.seq))
    return motifs

# Gets sequences from a fasta file. Calls findMotif for each.
def parse_fasta(fasta):
    all_motifs = []
    for seq_record in SeqIO.parse(fasta, "fasta"): #for each entry do the below.
        terminal_width = os.get_terminal_size()[0]
        print(Fore.LIGHTGREEN_EX)
        for x in range(terminal_width):
            print(Fore.LIGHTGREEN_EX + '-', end = '')
            
        print('\n')
        motifs = slice_motifs(str(seq_record.seq), seq_record.id)
        all_motifs.append(motifs)
        print_motifs(motifs)
    return all_motifs


def get_d1d1(start, end, seq):
    d1d1_search_area = seq[0:150] # Limits the area in which the d1d1 region is found, usually.
    d1d1_start_position = d1d1_search_area.find(start, 0, 20)# Find where the starting pattern is at. Limit to the first 20 bases.
    if (d1d1_start_position == -1):#If the start position is not found, return None to the main program
        return None
    else: # But if it is found, try searching for the end pattern.
        end_matches = re.finditer(end, d1d1_search_area)#Find all the matching bases to the pattern in the argument passed to the function (end)
        d1d1_results = []
        
        for match in end_matches: #For each match, add the end position to the d1d1_results array as a (start,end) tuple.
            d1d1_results.append(seq[d1d1_start_position:match.end()])
        return d1d1_results
    
#Sometimes the regex tends to find the first (longest) start match, but the real BoxB is shorter. There is a matching start pattern down the line, within the pattern match
#So this will cut the string by the start of the match + 1 so it doesn't find it again and looks for the next
def get_motif(start, end, seq, min_bases = 1, max_bases = 200):
    result = [] # Array where matches are stored
    pattern = r"(" + start + "(.*?)" + end + ")"
    while len(seq): #while the length of the sequence evaluated is not 0
        match = re.search(pattern, seq) #find the pattern
        if (match): #if found
            if(len(match.group()) > 15 ) & (len(match.group()) < 80 ): #if the length of the seq is <80
                result.append(str(match.group())) #add it to the result array

            seq = seq[match.start() + 1:] # slices the match off the seq for the next search.
        else:
            break 
    if len(result) == 0:
        return None
    return result



def slice_motifs(seq_input, organism_name): # returns a dictionary where each key is a motif and the value is the possible motifs.
    
    minimum_its_length = 20 #Used to filter out bad results.
    motifs = { organism_name : [],
                "leader" : [],
                "d1d1" : [],
                "sp_d2d3_sp" : [],
                "tRNA1" : [],
                "v2" : [],
                "tRNA2" : [],
                "BoxB" : [],
                "BoxA" : [],
                "D4" : [],
                "V3" : []
                }
    
    pico_cyano_flag = 0 #Used to identify picocyano d1d1 starts.

    its_seq_search = re.search(r"CCTCCTT", seq_input) #Extracting only the ITS region from a 16S-23S region. Find CCTCCTT, set ITS start position after that.
    if (its_seq_search == None):
        its_seq_search = re.search(r"CCTCCTA", seq_input)
        pico_cyano_flag = 1
        
    if ((its_seq_search == None)): 
        menu_option = ''
        print(Fore.RED + "Could not find the end of 16S to determine the ITS region boundaries. Results may be inaccurate.")
        its_start_position = 0
    else: #If the end of 16S was found, check if the length is too short.
        if (len(seq_input[its_seq_search.start():]) < minimum_its_length):
            print (Back.RED + Fore.WHITE + "Region length too short. Skipping this sequence.")
            return None
        its_start_position = its_seq_search.start() + 7 #Set the start position of the 
    
    its_seq = seq_input[its_start_position:]   #Sequence to be used for the motif search
    
    motifs[organism_name].append(its_seq) #store the whole ITS region to be used as a reference.
     
    if (pico_cyano_flag == 1): # Change D1D1 start if we are dealing with a picocyano. 
        d1d1 = get_d1d1("GACAA", r"[AT]TGTC", its_seq)
    else:
        d1d1 = get_d1d1("GACCT", r"AGGTC", its_seq)
        if (d1d1 == None):
            d1d1 = get_d1d1("GACCA", r"[AT]GGTC", its_seq)
        if (d1d1 == None): 
            d1d1 = get_d1d1("GACCG", r"[AC]GGTC", its_seq)
        if (d1d1 == None): 
            d1d1 = get_d1d1("GACCC", r"[AC]GGTC", its_seq)
            
    if (d1d1 == None):
            motifs["leader"] = None
            motifs["d1d1"] = None
    else: 
        leader_start = its_seq.find(d1d1[0]) #get the starting index of the d1d1
        motifs["leader"].append(its_seq[0:leader_start]) #use the above to define where leader ends (start of d1d1)
        for seq in d1d1: 
            motifs["d1d1"].append(seq) #Append the d1d1 results to the d1d1 dict key (an array).
        its_seq = its_seq[its_seq.rindex(motifs["d1d1"][-1]):] #Trim the its_seq. rindex picks the last index of the found string. Equivalent to the end of the d1d1 string. use [-1] to pick the longest d1d1 result (the last one found, latest index of the d1d1 list)
    
    
    #spacer-D2-spacer D3-spacer region and  tRNA-1
    trna1 = get_motif("GGGC[TC]ATTA","GGCCCA", its_seq)
    #trna1_results = re.search(r"GGGC[TC]ATTA(.*?)GGCCCA",its_seq) #find text between basal clamps
    if (trna1 == None):
        motifs["tRNA1"] = None
        motifs["sp_d2d3_sp"] = None
        #TODO: what to do with sp_d2d3_sp if tRNA1 is not found.
    else:  
        motifs["sp_d2d3_sp"].append(its_seq[0:its_seq.find(trna1[0])]) #Store what is between the end of d1d1 (first char) and the start of tRNA1
        for seq in trna1: 
            motifs["tRNA1"].append(seq)
        #motifs["tRNA1"] = trna1 #store tRNA1 sequence string found in its_seq                 
        its_seq = its_seq[its_seq.rindex(motifs["tRNA1"][-1]):] #trim processed its region
    
    #find tRNA-Ala(2) 
    trna2 = get_motif("GGGG", "[TC]CTCCA", its_seq)
    if (trna2 == None):
        motifs["tRNA2"] = None
        motifs["v2"] = None
    else: 
        motifs["v2"].append(its_seq[0:its_seq.find(trna2[0])]) #everything between end of tRNA1 and start of tRNA2
       # motifs["tRNA2"] = trna2
        for seq in trna2: 
            motifs["tRNA2"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA2"][-1]):] #trim processed its region
        
    #find BoxB
    boxb = get_motif("CAGC","GCTG", its_seq)
    if (boxb == None):
        boxb = get_motif("AGCA","CTG", its_seq)
    if (boxb == None):
        motifs["BoxB"] = None
    else:        
        for seq in boxb: 
            motifs["BoxB"].append(seq)
        #motifs["BoxB"] = boxb        
        #rindex picks the last index of the found string. The first BoxB [0] is the longest, so using that.
        its_seq = its_seq[its_seq.rindex(motifs["BoxB"][0]):]
        
    #find Box-A
    boxa = get_motif("TCAAA[GA]G","A[AC]G", its_seq[::-1])
    if (boxa == None):
        motifs["BoxA"] = None
    else: 
        for seq in boxa:
            motifs["BoxA"].append(seq[::-1])
        #motifs["BoxA"] = boxa[0][::-1] #[0] is the string of the result, [::-1] reverses it.
        its_seq = its_seq[its_seq.rindex(motifs["BoxA"][0]) - 3:]
        #boxa_end = re.search(motifs["BoxA"], its_seq) #Search again for string to find the end position. Since the string was reversed, the end() of the result above can't be used, search for the found string to find the real end position
        #its_seq = its_seq[boxa_end.end() - 3:] #trim remaining its region, leave ACT in the its_seq   
     
    #find D4  
    d4 = get_motif("ACT", "TA[TGAC]", its_seq)
    if (d4 == None):
        motifs["D4"] = None
    else:  
        for seq in d4: 
           motifs["D4"].append(seq)
        #motifs["D4"] = d4
        its_seq = its_seq[its_seq.rindex(motifs["D4"][0]):]
        
    #Find v3
    v3 = get_motif("GT[CA]G", "CA[CG]A[GC]", its_seq)
    if (v3 == None):
        motifs["V3"] = None
    else:  
        for seq in v3: 
           motifs["V3"].append(seq)
        #motifs["v3"] = v3
        its_seq = its_seq[its_seq.rindex(motifs["V3"][0]):]
        
    return motifs
   
# If -g then run fetch then parse_fasta, 
# If -f then just run parse_fasta

# Create the parser    
parser = argparse.ArgumentParser(description = "Find motifs within Cyanobacteria ITS region sequences")

group = parser.add_mutually_exclusive_group() #Create a group of mutually exclusive argument options

group.add_argument('-f', 
                   '--fasta', 
                   action = "store", 
                   help = "Find motifs in a given fasta file.",
                   ) #Add arguments to group
group.add_argument('-g', 
                   '--genbank', 
                   help = "Fetch a sequence from a given genbank accession number.",
                   nargs='+',
                   )
parser.add_argument('-m',
                    '--motif',
                    help = "Select which motifs to extract",
                    default = "all",
                    nargs="*", #Expects 0 or more values, if none, then default applies.
                    choices=('leader', 'd1d1', 'v2', 'trna1', 'trna2', 'boxa', 'boxb', 'd4', 'v3', 'all')
                    )
parser.add_argument('-e',
                    '--email',
                    help = "Provide email for genbank query. If not provided, you will be prompted for one anyway.",
                    default = "none",
)

#Process the arguments
if len(sys.argv) == 1: #If there is no argument then print the help
    parser.print_help()
    parser.exit()
args = parser.parse_args()#Parse the arguments, store them in args.

def check_email(email):
    valid_email_format = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
    if(re.fullmatch(valid_email_format, email)):
        return True
    else:
        return False
    
if(args.fasta):
    try:
        fasta = open(args.fasta, "r")
        print("Processing fasta file...")
        motifs = parse_fasta(fasta)
        generate_json(motifs)
        exit()
    except IOError:
        print ("File not found.")
        exit() 
    
if(args.genbank):
    if (args.email): #Get email to send to Entrez. Either from argument or ask the user for input
        email = args.email
    else:
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
    while (not check_email(email)):#Check that the email is valid. 
        print("Invalid email, try again: ")
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
        
    for gb in args.genbank:
        print("Fetching genbank data from " + gb)
        try:
            motifs = parse_genbank(gb, email)
            print_motifs(motifs)
        except:
            print("Error parsing that accession number. Exiting.")
            exit()