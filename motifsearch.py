from asyncore import file_dispatcher
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
                    print(Fore.CYAN + Style.BRIGHT + key + "\n" + 
                          Fore.MAGENTA +  "Sequence: " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(motifs[key]) + 
                          Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(motifs[key])) + "\n")    
# Gets a sequence from genbank
def parse_genbank(accession, email): #fetch sequence and taxonomy with accession#
    Entrez.email = email # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
        seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        motifs = findMotifs(str(seq_record.seq))
    return motifs
# Gets sequences from a fasta file. Calls findMotif for each.
def parse_fasta(fasta):
    for seq_record in SeqIO.parse(fasta, "fasta"): #for each entry do the below.
        print(Fore.LIGHTGREEN_EX + "\nOrganism name: " + seq_record.id + ":")
        print("\n")
        motifs = findMotifs(str(seq_record.seq))
        return motifs
    return

# Finds all the possible D1D1 sequences. Gets the start and end patterns from the main program.

#Function gets the start pattern, the end pattern, and the full sequence.
def get_d1d1(start, end, seq):
    
    # Limits the area in which the d1d1 region is found, usually.
    d1d1_search_area = seq[0:150]
    # Find where the starting pattern is at. Limit to the first 20 bases.
    d1d1_start_position = d1d1_search_area.find(start, 0, 20)
    
    #If the start position is not found, return None to the main program
    if (d1d1_start_position == -1):
        return None
    # But if it is found, try searching for the end pattern.
    else:
        #Find all the matching bases to the pattern in the argument passed to the function (end)
        end_matches = re.finditer(end, d1d1_search_area)
        d1d1_results = []
        #For each match, add the end position to the d1d1_results array as a (start,end) tuple.
        for match in end_matches:
            d1d1_results.append(seq[d1d1_start_position:match.end()])
            #d1d1_results.append((d1d1_start_position,match.end()))
        #return it back to the program.
        return d1d1_results
# returns a dictionary where each key is a motif and the value is the possible motifs.
def findMotifs(seq_input):
    
    minimum_its_length = 20 #Used to filter out short results.
    motifs = {} #dictionary. motifs[motif-name]=[start-position,end-position, sequence, length]
    pico_cyano_flag = 0 #Used to identify picocyano d1d1 starts.

    #Extracting only the ITS region from a 16S-23S region.
    # Find CCTCCTT, set ITS start position after that.
    its_seq_search = re.search(r"CCTCCTT", seq_input)
    if (its_seq_search == None):
        its_seq_search = re.search(r"CCTCCTA", seq_input)
        pico_cyano_flag = 1
        
    #If it cannot find the end of 16S, ask if they want to continue anyway in case they provided the ITS region only.
    if ((its_seq_search == None)): 
        menu_option = '' #Initializes the menu option.
        print(Fore.RED + "Could not find the end of 16S to determine the ITS region boundaries. Results may be inaccurate.")
        its_start_position = 0
    else: #If the end of 16S was found, check if the length is too short.
        if (len(seq_input[its_seq_search.start():]) < minimum_its_length):
            print (Back.RED + Fore.WHITE + "Region length too short. Skipping this sequence.")
            return None
        its_start_position = its_seq_search.start() + 7 #Set the start position of the 
    
    #Sequence to be used for the motif search
    its_seq = seq_input[its_start_position:] 
    
    # NAME : SEQUENCE
    motifs["ITS"] = its_seq #store the whole ITS region to be used as a reference.
     
    # Change D1D1 start if we are dealing with a picocyano. 
    if (pico_cyano_flag == 1):
        d1d1_results = get_d1d1("GACAA", r"[AT]TGTC", its_seq)
    else:
        d1d1_results = get_d1d1("GACCT", r"AGGTC", its_seq)
        if (d1d1_results == None):
            d1d1_results = get_d1d1("GACCA", r"[AT]GGTC", its_seq)
        if (d1d1_results == None): 
            d1d1_results = get_d1d1("GACCG", r"[AC]GGTC", its_seq)
    
    if (d1d1_results == None):
            motifs["leader"] = None
            motifs["d1d1"] = None
    #get_d1d1 returns an array with as many sequences as it found.
    else: 
        #get the starting index of the d1d1
        leader_start = its_seq.find(d1d1_results[0])
        #use the above to define where leader ends (start of d1d1)
        motifs["leader"] = its_seq[0:leader_start] #The sequence string of the motif
        #create an empty list for the d1d1 results
        motifs["d1d1"] = []
        #Append the d1d1 results to the d1d1 key, as an array.
        for seq in d1d1_results:
            motifs["d1d1"].append(seq)
        
        #Trim the its_seq. rindex picks the last index of the found string. Equivalent to the end of d1d1. use [-1] to pick the 
        # longest d1d1 result (the last one found, latest index of the d1d1 list)
        its_seq = its_seq[its_seq.rindex(motifs["d1d1"][-1]):] #trim processed its region to avoid overlapping searches.
    
    
    # spacer-D2-spacer D3-spacer regionand  tRNA-Ile(1)
    tRNA1Search = re.search(r"GGGC[TC]ATTA(.*?)GGCCCA",its_seq) #find text between basal clamps
    if (tRNA1Search == None):
        motifs["tRNA1"] = None
        #TODO: what to do with sp_d2d3_sp if tRNA1 is not found.
    else:  
        # leader and d1d1 were trimmed from its_seq above. The beginning of the remaining string is the start of spacer.
        motifs["sp_d2d3_sp"] = its_seq[0:tRNA1Search.start()] #Store what is between the end of d1d1 (first char) and the start of tRNA1
        motifs["tRNA1"] = tRNA1Search[0] #store tRNA1 sequence string found in its_seq                 
        its_seq = its_seq[tRNA1Search.end():] #trim processed its region
    
    #find tRNA-Ala(2)
    tRNA2Search = re.search(r"GGGG(.*?)[TC]CTCCA",its_seq)  #find text between basal clamps
    if (tRNA2Search == None):
        motifs["tRNA2"] = None
    else: 
        motifs["v2"] = its_seq[0:tRNA2Search.start()] #everything between end of tRNA1 and start of tRNA2
        motifs["tRNA2"] = tRNA2Search[0]
        its_seq = its_seq[tRNA2Search.end():] #trim processed its region
        
    #find BoxB
    if (pico_cyano_flag == 1):
        BoxBSearch = re.search(r"CAGC(.*?)GCTG", its_seq) #search for this is if's a pico
    else:
        BoxBSearch = re.search(r"[TC]AGCA[ACT](.*?)TGCT[AG]", its_seq) #find text between basal clamps
    if (BoxBSearch == None):
        motifs["BoxB"] = None
    else:        
        
        motifs["BoxB"] = BoxBSearch[0]
        its_seq = its_seq[BoxBSearch.end():] #trim remaining its region
        
    #find Box-A
    BoxASearch = re.search(r"TCAAA[GA]G(.*?)A[AC]G", its_seq[::-1]) #BoxASearch[0] includes the string reversed.
    if (BoxASearch == None):
        motifs["BoxA"] = None
    else: 
        motifs["BoxA"] = BoxASearch[0][::-1] #[0] is the string of the result, [::-1] reverses it.
        #Search again for string to find the end position.
         #Since the string was reversed, the end() of the result above can't be used, search for the found string to find the real end position
        BoxAEndSearch = re.search(motifs["BoxA"], its_seq)
        its_seq = its_seq[BoxAEndSearch.end() - 3:] #trim remaining its region, leave ACT in the its_seq   
     
    #find D4  
    D4Search = re.search(r"ACT(.*?)TA[TGAC]",its_seq) #Should be right after the GAA above
    if (D4Search == None):
        motifs["D4"] = None
    else:  
        motifs["D4"] = D4Search[0]
        its_seq = its_seq[D4Search.end():] #trim remaining its region
        
    #Find v3
    V3Search = re.search(r"GT[CA]G(.*?)CA[CG]A[GC]",its_seq) 
    if (V3Search == None):
        motifs["v3"] = None
    else:  
        motifs["v3"] = V3Search[0]
        its_seq = its_seq[V3Search.end():]
        
    return motifs
    # end of extractMotif
   
   
# If -g then run fetch then parse_fasta, 
# If -f then just run parse_fasta

# Create the parser    
parser = argparse.ArgumentParser(description = "Find motifs within Cyanobacteria ITS region sequences")

#Create a group of mutually exclusive argument options
group = parser.add_mutually_exclusive_group()

#Add arguments to group
group.add_argument('-f', 
                   '--fasta', 
                   action = "store", 
                   help = "Find motifs in a given fasta file.",
                #    type = argparse.FileType('r'),
                   )
group.add_argument('-g', 
                   '--genbank', 
                   help = "Fetch a sequence from a given genbank accession number.",
                   nargs='+',
                   )
# Adds standalone arguments
# parser.add_argument('-H', 
#                     '--html',  
#                     help = "Outputs an HTML file with motifs highlighted over the ITS sequences. Optional: supply a filename after the flag. Default: its_output.html",
#                     default = "its_output.html",
#                     type = "str",
#                     )
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

#If there is no argument then print the help
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

#Parse the arguments, store them in args.
args = parser.parse_args()

def check_email(email):
    valid_email_format = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
    if(re.fullmatch(valid_email_format, email)):
        return True
    else:
        return False
    
if(args.fasta):
    #TODO: make sure the file exists.
    try:
        fasta = open(args.fasta, "r")
        print("Processing fasta file...")
        motifs = parse_fasta(fasta)
        print_motifs(motifs)
        generate_json(motifs)
        exit()
    except IOError:
        print ("File not found.")
        exit() 
    
if(args.genbank):
    #Get email to send to Entrez. Either from argument or ask the user for input
    if (args.email):
        email = args.email
    else:
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
    #Check that the email is valid. 
    while (not check_email(email)):
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

# if(args.html):
#     generate_html(args.html, args.motif)
    
    

    