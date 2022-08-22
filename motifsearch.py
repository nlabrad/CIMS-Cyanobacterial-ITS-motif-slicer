import os
import re
import sys
import argparse
import colorama
import json
from Bio import Entrez, SeqIO #import the required libraries
from colorama import Fore, Back, Style
colorama.init(autoreset=True)

#TODO: add logging.
#TODO: option to create CSV
#TODO: HTML with highlighted motif text
#TODO: option to select just one motif from given seq instead of all the motifs.
#TODO: option to exclude ones without tRNAs, or alternatively, just include those with both tRNAs

def generate_json(motif_list):
    """Generates json file from the dictionary argument (usually the motifs)

    Args:
        motif_list (dictionary): motifs
    """
    
    json_string = json.dumps(motif_list, indent=3)
    file = open("motifs.json", "w", encoding='utf-8')
    file.write(json_string)
    file.close()
    
def check_email(user_email):
    """ Validates email format

    Args:
        email (string): email to validate

    Returns:
        boolean: True if valid, False if invalid.
    """
    valid_email_format = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
    return re.fullmatch(valid_email_format, user_email)
    
def print_motifs(motif_list):
    """Gets the dict of motifs and prints them out in the terminal

    Args:
        motifs (dict): list of motifs to print
    """
    if motif_list is None:
        print(Back.RED + Fore.WHITE + "ITS Region not found in this sequence!")
    else:
        for key in motif_list:
            if motif_list[key] is None:
                if (key == "tRNA1" or key == "tRNA2"):
                    print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not present in this operon.")
                    print("\n")  
                else: 
                    print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not found in this sequence.")
                    print("\n")
            else:
                if len(motif_list[key]) > 1:
                    print(Fore.CYAN + Style.BRIGHT + key + " " +Fore.RED + Style.BRIGHT + str(len(motif_list[key])) + " possible sequences found: \n")     
                    for index, item in enumerate(motif_list[key]):
                        print(Fore.MAGENTA + Style.BRIGHT +  "Sequence " + str(index+1) + ": " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                            Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")
                else:
                    print(Fore.CYAN + Style.BRIGHT + key + ":")
                    for item in motif_list[key]:
                        print(Fore.MAGENTA + Style.BRIGHT +  "Sequence " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                            Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")

# Gets the specified sequence from genbank to process.
def parse_genbank(accession, valid_email): #fetch sequence and taxonomy with accession#
    """ Gets fasta file from genbank to be processed

    Args:
        accession (string): accession number to get from genbank
        email (string): email address to be used when contacting GB

    Returns:
        dict: list of motifs
    """
    Entrez.email = valid_email # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
        seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        sliced_motifs = slice_motifs(str(seq_record.seq), str(seq_record.id))
    return sliced_motifs

# Gets sequences from a fasta file. Calls findMotif for each.
def parse_fasta(fasta_string):
    """ Processes a given fasta file.

    Args:
        fasta (string): fasta file as a string

    Returns:
        list: list with all the organisms in the fasta file
    """
    all_motifs = []
    for seq_record in SeqIO.parse(fasta_string, "fasta"): #for each entry do the below.
        terminal_width = os.get_terminal_size()[0]
        print(Fore.LIGHTGREEN_EX)
        for x in range(terminal_width):
            print(Fore.LIGHTGREEN_EX + '-', end = '')
            
        print('\n')
        sliced_motifs = slice_motifs(str(seq_record.seq), seq_record.id)
        all_motifs.append(sliced_motifs)
        print_motifs(sliced_motifs)
    return all_motifs


def get_d1d1(start, end, seq):
    """ Algorithm to search for d1d1 motif

    Args:
        start (string): d1d1 start pattern
        end (string): d1d1 end pattern
        seq (string): sequence where the d1d1 is to be found.

    Returns:
        list: list containing all the possible d1d1 sequences.
    """
    d1d1_search_area = seq[0:150] # Limits the area in which the d1d1 region is found, usually.
    d1d1_start_position = d1d1_search_area.find(start, 0, 20)# Find where the starting pattern is at. Limit to the first 20 bases.
    if d1d1_start_position == -1: #If the start position is not found, return None to the main program
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
    
    """ Function that searches for a motif given the parameters
    
    Args:
        start (string): motif start pattern
        end (string): motif end pattern
        seq (string): sequence where the d1d1 is to be found.
        min_bases (int): OPTIONAL argument to set the minimum length of results
        max_bases (int): OPTIONAL argument to set the maximum length of results

    Returns:
        list: list of sequences (string) that correspond to the motif
    """
    result = [] # Array where matches are stored
    pattern = r"(" + start + "(.*?)" + end + ")"
    while len(seq): #while the length of the sequence evaluated is not 0
        match = re.search(pattern, seq) #find the pattern
        if match: #if found
            if(len(match.group()) > min_bases ) & (len(match.group()) < max_bases): #if the length of the seq is <80
                result.append(str(match.group())) #add it to the result array

            seq = seq[match.start() + 1:] # slices the match off the seq for the next search.
        else:
            break 
    if len(result) == 0:
        return None
    return result

def slice_its_region(seq_input, min_length):
    """ Slices only the ITS region from the sequence given by the user (fasta or genbank)

    Args:
        seq_input (string): full sequence
        min_length (minimum length): establishes the minimum length that the ITS region can be to avoid false positives

    Returns:
        int: start position of the ITS region (index relative to the raw sequence provided)
    """
    its_seq_search = re.search(r"CCTCCTT", seq_input)
    
    if its_seq_search is None:
        its_seq_search = re.search(r"CCTCCTA", seq_input)
        global PICO_CYANO_FLAG
        PICO_CYANO_FLAG = 1
    if its_seq_search is None:
        print(Fore.RED + "Could not find the end of 16S to determine the ITS region boundaries. Results may be inaccurate.")
        return 0
    else:
        if len(seq_input[its_seq_search.start():]) < min_length:
            print (Back.RED + Fore.WHITE + "Region length too short. Skipping this sequence.")
            return None
        return (its_seq_search.start() + 7)
    
def slice_motifs(seq_input, organism_name):
    """Main function that coordinates the calls to find all the motifs. Contains the motif patterns

    Args:
        seq_input (string): sequence where to search for the motifs
        organism_name (string): organism name pulled from the first line of the fasta file.

    Returns:
        dict: dictionary with motifs. Key: motif name, value: motif sequence(s) (list)
    """
    minimum_its_length = 20#Used to filter out bad results.
    its_start_position = 0#ITS region start index relative to the fasta/genbank supplied sequence
    global PICO_CYANO_FLAG #Used to identify picocyano d1d1 starts.
    PICO_CYANO_FLAG = 0
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
    # its_seq_search = re.search(r"CCTCCTT", seq_input)
    # if its_seq_search is None:
    #     its_seq_search = re.search(r"CCTCCTA", seq_input)
    #     pico_cyano_flag = 1
        
    # if its_seq_search is None:
    #     print(Fore.RED + "Could not find the end of 16S to determine the ITS region boundaries. Results may be inaccurate.")
        
    # else:
    #     if len(seq_input[its_seq_search.start():]) < minimum_its_length:
    #         print (Back.RED + Fore.WHITE + "Region length too short. Skipping this sequence.")
    #         return None
    #     its_start_position = its_seq_search.start() + 7
    its_start_position = slice_its_region(seq_input, minimum_its_length)
    if its_start_position is None:
        print (Back.RED + Fore.WHITE + "Found ITS Region length is too short (not accurate). Skipping this sequence.")
        return None
        
    its_seq = seq_input[its_start_position:]   #Sequence to be used for the motif search. Found motifs get removed from the seq before moving on to the next one.
    motifs[organism_name].append(its_seq) #store the whole ITS region sequence in the dict. This one does not change with the search.
    
    if PICO_CYANO_FLAG == 1: # Change D1D1 start if we are dealing with a picocyano. 
        d1d1 = get_d1d1("GACAA", r"[AT]TGTC", its_seq)
    else:
        d1d1 = get_d1d1("GACCT", r"AGGTC", its_seq)
        if d1d1 is None:
            d1d1 = get_d1d1("GACCA", r"[AT]GGTC", its_seq)
        if d1d1 is None:
            d1d1 = get_d1d1("GACCG", r"[AC]GGTC", its_seq)
        if d1d1 is None:
            d1d1 = get_d1d1("GACCC", r"[AC]GGTC", its_seq)
            
    if d1d1 is None:
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
    if trna1 is None:
        motifs["tRNA1"] = None
        motifs["sp_d2d3_sp"] = None
    else:
        motifs["sp_d2d3_sp"].append(its_seq[0:its_seq.find(trna1[0])]) #Store what is between the end of d1d1 (first char) and the start of tRNA1
        for seq in trna1:
            motifs["tRNA1"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA1"][-1]):] #trim processed its region
    
    #find tRNA-Ala(2)
    trna2 = get_motif("GGGG", "[TC]CTCCA", its_seq)
    if trna2 is None:
        motifs["tRNA2"] = None
        motifs["v2"] = None
    else: 
        motifs["v2"].append(its_seq[0:its_seq.find(trna2[0])]) #everything between end of tRNA1 and start of tRNA2
        for seq in trna2:
            motifs["tRNA2"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA2"][-1]):] #trim processed its region
        
    #find BoxB
    boxb = get_motif("CAGC","GCTG", its_seq, 1, 80)
    if boxb is None:
        boxb = get_motif("AGCA","CTG", its_seq)
    if boxb is None:
        motifs["BoxB"] = None
    else:        
        for seq in boxb:
            motifs["BoxB"].append(seq)
        #rindex picks the last index of the found string. The first BoxB [0] is the longest, so using that.
        its_seq = its_seq[its_seq.rindex(motifs["BoxB"][0]):]
        
    #find Box-A
    boxa = get_motif("TCAAA[GA]G","A[AC]G", its_seq[::-1], 1, 80) #reverse the string, look from the back to the front
    if boxa is None:
        motifs["BoxA"] = None
    else:
        for seq in boxa:
            motifs["BoxA"].append(seq[::-1]) #reverses it back.
        its_seq = its_seq[its_seq.rindex(motifs["BoxA"][0]) - 3:]
        
    #find D4
    d4 = get_motif("ACT", "TA[TGAC]", its_seq)
    if d4 is None:
        motifs["D4"] = None
    else:
        for seq in d4:
            motifs["D4"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["D4"][0]):]
        
    #Find v3
    v3 = get_motif("GT[CA]G", "CA[CG]A[GC]", its_seq)
    if v3 is None:
        motifs["V3"] = None
    else:
        for seq in v3:
            motifs["V3"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["V3"][0]):]
        
    return motifs


#Main
# Create the parser
parser = argparse.ArgumentParser(
    description = "Find motifs within Cyanobacteria ITS region sequences")

group = parser.add_mutually_exclusive_group()#Create a group of mutually exclusive argument options

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
                    help = "Provide email for genbank query. If not provided, you will be prompted for one at runtime.",
                    default = None,
)

#Process the arguments
if len(sys.argv) == 1: #If there is no argument then print the help
    parser.print_help()
    parser.exit()
args = parser.parse_args()#Parse the arguments, store them in args.


if args.fasta:
    try:
        fasta = open(args.fasta, "r", encoding="UTF-8")
        print("Processing fasta file...")
        motifs = parse_fasta(fasta)
        generate_json(motifs)
        exit()
    except IOError:
        print ("File not found.")
        exit()
    
if args.genbank:
    if args.email: #Get email to send to Entrez. Either from argument or ask the user for input
        email = args.email
    else:
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
    while (not check_email(email)): #Check that the email is valid.
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
            