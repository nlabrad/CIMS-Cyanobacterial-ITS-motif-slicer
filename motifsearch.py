import re, sys, argparse, colorama
from Bio import Entrez, SeqIO #import the required libraries
from colorama import Fore, Back, Style
colorama.init(autoreset=True)

#TODO: add logging.
#TODO: handle not having an argument after the option -f or -g
#TODO: handle not finding the file passed.

# Gets a sequence from genbank
def parse_genbank(accession, email): #fetch sequence and taxonomy with accession#
    Entrez.email = email # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype="fasta", retmode="text") as handle: # id is what genbank ID to query, type is genbank
        seq_record = SeqIO.read(handle, "fasta") #get sequence in fasta format to be used as needed.
        motifs = findMotifs(str(seq_record.seq))
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
                    print(Fore.CYAN + Style.BRIGHT + key + "\n" + Fore.MAGENTA +  "Sequence: " + Fore.LIGHTYELLOW_EX + Style.NORMAL + motifs[key][2] + Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(motifs[key][3]) + "\n")
    return


# Gets sequences from a fasta file. Calls findMotif for each.
def parse_fasta(fasta):
    for seq_record in SeqIO.parse(fasta, "fasta"): #for each entry do the below.
        print(Fore.LIGHTGREEN_EX + "\nOrganism name: " + seq_record.id + ":")
        print("\n")
        motifs = findMotifs(str(seq_record.seq))
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
                    print(Fore.CYAN + Style.BRIGHT + key + "\n" + Fore.MAGENTA +  "Sequence: " + Fore.LIGHTYELLOW_EX + Style.NORMAL + motifs[key][2] + Style.BRIGHT + Fore.LIGHTMAGENTA_EX + " \nLength: " + Style.NORMAL + Fore.LIGHTYELLOW_EX + str(motifs[key][3]) + "\n")
    return
    
def findMotifs(seq): #Find the motifs
    motifs = {} #dictionary. motifs[motif-name]=[start-position,end-position, sequence, length]
    index_shift = 0

    #Extraction of ITS region from a 16S-23S region.
    # Find CCTCCTT, get rid of what's before that.
    its_seq_search = re.search(r"CCTCCTT", seq)
    
    if ((its_seq_search == None)):
        ans = ''
        while not (ans == 'N' or ans == 'Y'):
            print(Fore.RED + "Could not find the end of 16S to determine the ITS region boundaries.")
            ans = input(Fore.RED + "Proceed with search anyway? (Y/N): ").upper()
            if(ans == 'Y'):
                print("Proceeding with the whole sequence...\n")
                itsStartPosition = 0
                break
            if (ans == 'N'):
                print("Skipping this organism.\n")
                return None
            else:
                print(Fore.RED + "Invalid option. Valid Options: Y or N\n")
    else:
        if (len(seq[its_seq_search.start():]) < 20 and its_seq_search is not None):
            print (Back.RED + Fore.WHITE + "Region length too short. Skipped.")
            return None
        itsStartPosition = re.search("CCTCCTT", seq).start() + 7 #If it's >20
    
    its_region = seq[itsStartPosition:]
    
    motifs["ITS"] = [0, len(its_region)-1, its_region, len(its_region)] #store the whole ITS region to be used as a reference.
     
    d1d1Search = re.search(r"GACCT(.*?)AGGTC", its_region) #find text between basal clamps, starting with GACCT/C to the first AGGTC (*? is lazy search)
    if (d1d1Search == None): 
        d1d1Search = re.search(r"GACCA(.*?)[AT]GGTC",its_region) 
    if (d1d1Search == None): 
        d1d1Search = re.search(r"GACCG(.*?)CGGTC",its_region)
    if (d1d1Search == None):
            motifs["leader"] = None
            motifs["d1d1"] = None
    else: 
        motifs["leader"] = [0,d1d1Search.start()] #start-end positions relative to the ITS region
        motifs["leader"].append(motifs["ITS"][2][motifs["leader"][0]:motifs["leader"][1]]) #The sequence string of the motif
        motifs["leader"].append(len(motifs["leader"][2])) #The length of the motif
        
        motifs["d1d1"] = [d1d1Search.start(),d1d1Search.end()] #[0] = start, [1] = end
        motifs["d1d1"].append(motifs["ITS"][2][motifs["d1d1"][0]:motifs["d1d1"][1]]) #[2] = sequence string
        motifs["d1d1"].append(len(motifs["d1d1"][2])) #[3] = lenght
        
        index_shift = d1d1Search.end() #count of processed characters. Used to keep track of absolute position in the index.
        its_region = its_region[d1d1Search.end():] #trim processed its region to avoid overlapping searches.
    
    
    # spacer-D2-spacer D3-spacer regionand  tRNA-Ile(1)
    tRNA1Search = re.search(r"GGGC[TC]ATTA(.*?)GGCCCA",its_region) #find text between basal clamps
    if (tRNA1Search == None):
        motifs["tRNA1"] = None
    else:  
        motifs["sp_d2d3_sp"] = [d1d1Search.end(), tRNA1Search.start()+index_shift]
        motifs["sp_d2d3_sp"].append(motifs["ITS"][2][motifs["sp_d2d3_sp"][0]:motifs["sp_d2d3_sp"][1]])
        motifs["sp_d2d3_sp"].append(len(motifs["sp_d2d3_sp"][2]))
        
        motifs["tRNA1"] = [tRNA1Search.start()+index_shift, tRNA1Search.end()+index_shift]
        motifs["tRNA1"].append(motifs["ITS"][2][motifs["tRNA1"][0]:motifs["tRNA1"][1]])
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
        motifs["tRNA2"].append(motifs["ITS"][2][motifs["tRNA2"][0]:motifs["tRNA2"][1]])
        motifs["tRNA2"].append(len(motifs["tRNA2"][2])) 
        
        index_shift = index_shift + tRNA2Search.end()
        its_region = its_region[tRNA2Search.end():] #trim processed its region
        
    #find BoxB
    BoxBSearch = re.search(r"[TC]AGCA[AC](.*?)TGCT[AG]", its_region) #find text between basal clamps
    if (BoxBSearch == None):
        motifs["BoxB"] = None
    else:        
        motifs["BoxB"] = [BoxBSearch.start()+index_shift, BoxBSearch.end()+index_shift]
        motifs["BoxB"].append(motifs["ITS"][2][motifs["BoxB"][0]:motifs["BoxB"][1]])
        motifs["BoxB"].append(len(motifs["BoxB"][2]))
        index_shift = index_shift + BoxBSearch.end()
        its_region = its_region[BoxBSearch.end():] #trim remaining its region
        
    #find Box-A
    #Old version
    #BoxASearch = re.search(r"G[CA]A(.*?)G[AG]AAACT",its_region) #TODO: make it search backwards from ACT and find the first pattern matched.
    BoxASearch = re.search(r"TCAAA[GA]G(.*?)A[AC]G", its_region[::-1]) #BoxASearch[0] includes the string reversed.
    if (BoxASearch == None):
        motifs["BoxA"] = None
    else: 
        motifs["BoxA"] = [BoxASearch.start()+index_shift, BoxASearch.end()- 3 +index_shift] #-3 to exclude ACT
        motifs["BoxA"].append(BoxASearch[0][::-1]) #[0] is the string of the result, [::-1] reverses it.
        motifs["BoxA"].append(len(motifs["BoxA"][2]))
        index_shift = index_shift + BoxASearch.end() - 3
        its_region = its_region[BoxASearch.end() - 3:] #trim remaining its region   
    
    
    #find D4  
    D4Search = re.search(r"ACT(.*?)TA[TGAC]",its_region) #Should be right after the GAA above
    if (D4Search == None):
        motifs["D4"] = None
    else:  
        motifs["D4"] = [D4Search.start()+index_shift, D4Search.end()+index_shift]
        motifs["D4"].append(motifs["ITS"][2][motifs["D4"][0]:motifs["D4"][1]])
        motifs["D4"].append(len(motifs["D4"][2]))
        index_shift = index_shift + D4Search.end()
        its_region = its_region[D4Search.end():] #trim remaining its region
        
    #Find v3
    V3Search = re.search(r"GT[CA]G(.*?)CA[CG]A[GC]",its_region) 
    if (V3Search == None):
        motifs["v3"] = None
    else:  
        motifs["v3"] = [V3Search.start()+index_shift, V3Search.end()+index_shift]
        motifs["v3"].append(motifs["ITS"][2][motifs["v3"][0]:motifs["v3"][1]])
        motifs["v3"].append(len(motifs["v3"][2]))
        index_shift = index_shift + V3Search.end()
        its_region = its_region[V3Search.end():]
        
    return motifs
    # end of extractMotif
   
def generate_html(file, motifs):
    f = open(file, 'w')
    html_header = """
    <html>
        <head>
        <title>Motif Search Results</title>
        <head>
        
        <body>
        <h1> Motif Search Results </h1>"""
    
    html_footer = """
    </body>
    </html>"""
    
    f.write(html_header)
       
    html_body = """
    <p>Test</p>
    """
    f.write(html_body)
    f.write(html_footer)
    f.close()
   
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
                   help = "Fetch a sequence from a given genbank accession number."
                   )

#Add standalone argument
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
        parse_fasta(fasta)
        exit()
    except IOError:
        print ("File not found")
        exit() 
    
if(args.genbank):
    email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
    valid_email = check_email(email)
    while (not valid_email):
        print("Invalid email, try again: ")
        email = input("Entrez requires an email to be associated with the requests. Please enter your email: ".lower())
        valid_email = check_email(email)
    
    print("Fetching genbank data...")
    try:
        parse_genbank(args.genbank, email)
    except:
        print("Error parsing that accession number. Exiting.")
        exit()

# if(args.html):
#     generate_html(args.html, args.motif)
    
    

    