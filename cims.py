import os
import re
import sys
import argparse
import json
import time
import colorama

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

def print_motifs(motif_list, print_list):
    """Gets the dict of motifs and prints them out in the terminal

    Args:
        motifs (dict): list of motifs to print
    """
    if motif_list is None:
        print(Back.RED + Fore.WHITE + "ITS Region not found in this sequence!")
    else:
        for key in motif_list:
            if (str(key) in print_list or "all" in print_list):
                if motif_list[key] is None: #If the key is empty
                    if (key == "tRNA_ile" or key == "tRNA_ala"):
                        print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not present in this operon.")
                        print("\n")
                    else:
                        print(Fore.LIGHTCYAN_EX + Style.BRIGHT + key + Fore.RED + " Not found in this sequence.")
                        print("\n")
                else:
                    if len(motif_list[key]) > 1: #Check if there is more than one sequence in that key
                        print(Fore.CYAN + Style.BRIGHT + key + " \n\t" +Fore.RED + Style.BRIGHT + str(len(motif_list[key])) + " possible sequences found! \n")     
                        for index, item in enumerate(motif_list[key]):
                            print(Fore.MAGENTA + Style.BRIGHT +  "\tSequence " + str(index+1) + ": " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                                Style.BRIGHT + Fore.LIGHTMAGENTA_EX + "\n\tLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")
                    else:
                        if key is list(motif_list.keys())[0]:
                            print(Fore.CYAN + Style.BRIGHT + key + " ITS Region:")
                        else:
                            print(Fore.CYAN + Style.BRIGHT + key + ":")
                        for item in motif_list[key]:
                            print(Fore.MAGENTA + Style.BRIGHT +  "\tSequence " + Fore.LIGHTYELLOW_EX + Style.NORMAL + str(item) + 
                                Style.BRIGHT + Fore.LIGHTMAGENTA_EX + "\n\tLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(len(item)) + "\n")

def parse_genbank(accession, valid_email): #fetch sequence and taxonomy with accession#
    """ Gets fasta file from genbank to be processed

    Args:
        accession (string): accession number to get from genbank
        email (string): email address to be used when contacting GB

    Returns:
        dict: list of motifs
    """
    Entrez.email = valid_email # email reported to entrez to associate with the query
    try:
        with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
            seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        sliced_motifs = slice_motifs(str(seq_record.seq), str(seq_record.id))
    except IOError:
        print("Network error. Could not fetch from genbank")
        raise IOError
    return sliced_motifs

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

def get_motif(start_pattern, end_pattern, seq, min_length=1, max_length=200):
    
    """Searches for the possible motifs given basal clamps in the arguments.
    Args: 
        start_pattern(string): motif opening basal clamp
        end_pattern(string): motif closing basal clamp
        seq(string): master sequence where to look for the motifs.
        min(int): minimum length of motif (default 1)
        max(int): maxmum length of motif (default 200)

    Returns:
        list: list of possible sequences.
    """
    motif_results = []
    while len(seq):
        start_match = re.search(start_pattern, seq)
        if start_match is None: #If the start position is not found, return None to the main program
            break
        seq = seq[start_match.start():]
        end_matches = re.finditer(r"" + end_pattern, seq)#Find all the matching bases to the pattern in the argument passed to the function (end)
        for match in end_matches: #For each match, add the end position to the d1d1_results array as a (start,end) tuple.
            motif_seq = str(seq[:match.end()])
            if len(motif_seq) > min_length and len(motif_seq) < max_length:
                motif_results.append(motif_seq)
        seq = seq[3:]
    if len(motif_results) == 0:
        return None
    return motif_results

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
        return its_seq_search.start() + 7

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
                "tRNA_ile" : [],
                "sp_v2_sp" : [],
                "tRNA_ala" : [],
                "BoxB" : [],
                "BoxA" : [],
                "D4" : [],
                "V3" : []
                }

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
    trna_ile = get_motif("GGGC[TC]ATTA","GGCCCA", its_seq, 70, 80)
    #trna1_results = re.search(r"GGGC[TC]ATTA(.*?)GGCCCA",its_seq) #find text between basal clamps
    if trna_ile is None:
        motifs["tRNA_ile"] = None
        motifs["sp_d2d3_sp"] = None
    else:
        motifs["sp_d2d3_sp"].append(its_seq[0:its_seq.find(trna_ile[0])]) #Store what is between the end of d1d1 (first char) and the start of tRNA_ile
        for seq in trna_ile:
            motifs["tRNA_ile"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA_ile"][-1]):] #trim processed its region

    #find tRNA-Ala(2)
    trna_ala = get_motif("GGGG", "[TC]CTCCA", its_seq, 70, 80)
    if trna_ala is None:
        motifs["tRNA_ala"] = None
        motifs["sp_v2_sp"] = None
    else: 
        motifs["sp_v2_sp"].append(its_seq[0:its_seq.find(trna_ala[0])]) #everything between end of tRNA_ile and start of tRNA_ala
        for seq in trna_ala:
            motifs["tRNA_ala"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA_ala"][-1]):] #trim processed its region

    #find BoxB
    boxb = get_motif("CAGC","GCTG", its_seq, 10, 50)
    if boxb is None:
        boxb = get_motif("AGCA","CTG", its_seq, 10, 50)
    if boxb is None:
        motifs["BoxB"] = None
    else:
        for seq in boxb:
            motifs["BoxB"].append(seq)
        #rindex picks the last index of the found string. The first BoxB [0] is the longest, so using that.
        its_seq = its_seq[its_seq.rindex(motifs["BoxB"][0]):]

    #find Box-A
    boxa = get_motif("TCAAA[GA]G","A[AC]G", its_seq[::-1], 1, 10) #reverse the string, look from the back to the front
    if boxa is None:
        motifs["BoxA"] = None
    else:
        for seq in boxa:
            motifs["BoxA"].append(seq[::-1]) #reverses it back.
        its_seq = its_seq[its_seq.rindex(motifs["BoxA"][0]) - 3:]

    #find D4
    d4 = get_motif("ACT", "TA[TGAC]", its_seq, 4, 20)
    if d4 is None:
        motifs["D4"] = None
    else:
        motifs["D4"].append(d4[0]) #D4 is not very accurate, store only the first one.
        its_seq = its_seq[its_seq.rindex(motifs["D4"][0]):]

    #Find v3
    v3 = get_motif("GT[CA]G", "CA[CG]A[GC]", its_seq, 4)
    if v3 is None:
        motifs["V3"] = None
    else:
        motifs["V3"].append(v3[0]) #V3 not very accurate, stores only the first one.
        its_seq = its_seq[its_seq.rindex(motifs["V3"][0]):]

    return motifs

global PICO_CYANO_FLAG
PICO_CYANO_FLAG = 0
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
parser.add_argument('-s',
                    '--select',
                    help = "Select which motifs to extract",
                    default = "all",
                    nargs="*", #Expects 0 or more values, if none, then default applies.
                    choices=('leader', 'd1d1', 'sp_v2_sp', 'trna1', 'trna2', 'boxa', 'boxb', 'd4', 'v3', 'all')
                    )
parser.add_argument('-e',
                    '--email',
                    help = "Provide email for genbank query. If not provided, you will be prompted for one at runtime.",
                    default = None,
                    )
parser.add_argument('-j',
                    '--json',
                    help = 'Outputs a json file with all the motifs',
                    action='store_true'
                    )
parser.add_argument('-o',
                   '--operon',
                   help = 'Returns how many tRNAs are found to use for homologous operon verification.',
                   action='store_true'
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
        if args.select:
            print_motifs(motifs, args.select)
        else:
            print_motifs(motifs, args.select)
        if args.json:
            generate_json(motifs) 
                
                    
    except IOError:
        print ("File not found.")
        exit()

if args.genbank:
    if args.email: #Get email to send to Entrez. Either from argument or ask the user for input
        email = args.email
    else:
        print("To query the Entrez database via a script, NCBI requires attaching an email address to the query.")
        print("This script does not check if the email is valid, does not store it and does not share it with anyone else.")
        print("The Biopython library is used to handle the communication with Entrez and the handling of the provided email address")
        email = input("Please enter your email for the Entrez query: ".lower())
    while (not check_email(email)): #Check that the email is valid.
        print("Invalid email, try again: ")
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
    for gb in args.genbank:
        try:
            print("\nFetching genbank data from " + gb + "\n")
            time.sleep(1) #Required to not go over the 3 queries/second threshold imposed by Entrez
            motifs = parse_genbank(gb, email)
            if args.select:
                print_motifs(motifs, args.select)
            else:
                print_motifs(motifs, args.select)
            if args.json:
                generate_json(motifs)
        except IOError:
            print("Error while parsing accession number. Exiting.")
            exit()