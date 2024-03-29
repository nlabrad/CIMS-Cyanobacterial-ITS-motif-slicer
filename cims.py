"""
Slices motifs out of an ITS Region sequence
"""

import argparse
from os import get_terminal_size
from pathlib import Path
from re import finditer, fullmatch, search
from sys import argv
from sys import exit as system_exit
from time import sleep
import colorama
from Bio import Entrez, SeqIO  # import the required libraries
from colorama import Back, Fore, Style

colorama.init(autoreset=True)


def parse_genbank(accession, email):
    """Gets fasta file from genbank to be processed

    Args:
        accession (string): accession number to get from genbank
        email (string): email address to be used when contacting GB

    Returns:
        dict: list of motifs
    """
    Entrez.email = email  # email reported to entrez to associate with the query
    try:
        # id is what genbank ID to query, type is genbank
        with Entrez.efetch(
            db="nucleotide", id=accession, rettype="fasta", retmode="text"
        ) as handle:
            seq_record = SeqIO.read(
                handle, "fasta"
            )  # get sequence in fasta format to be used as needed.
            handle.close()
        return seq_record
    except IOError:
        print ("Network error. Could not fetch from genbank")
    finally:
        system_exit(IOError)
    


def valid_email(unchecked_email):
    """Validates email format
    Args:
        email (string): email to validate

    Returns:
        boolean: True if valid, False if invalid.
    """
    valid_email_format = r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b"
    return fullmatch(valid_email_format, unchecked_email)


def parse_fasta(fasta_string):
    """Processes a given fasta file.

    Args:
        fasta (string): fasta file as a string

    Returns:
        list: list with all the organisms in the fasta file
    """
    seq_list = []
    for seq_record in SeqIO.parse(
        fasta_string, "fasta"
    ):  # for each entry do the below.
        seq_list.append(seq_record)
    return seq_list


def process_seqs(sequence_list):
    """Process the SeqRecords extracted from a fasta file or multiple genbank entries.

    Args:
        sequence_list (list): list of SeqRecords

    Returns:
        list: list of dictionaries that contain the motifs found in each SeqRecord
    """
    sliced_organisms = []
    for seq_record in sequence_list:  # for each entry do the below.
        # store the sliced motifs, pass the sequence and the organism name
        sliced_motifs = slice_motifs(str(seq_record.seq), seq_record.id)
        sliced_organisms.append(sliced_motifs)
    return sliced_organisms


def slice_motifs(seq_input, organism_name):
    """Main function that coordinates the calls to find all the motifs. Contains the motif patterns

    Args:
        seq_input (string): sequence where to search for the motifs
        organism_name (string): organism name pulled from the first line of the fasta file.

    Returns:
        dict: dictionary with motifs. Key: motif name, value: motif sequence(s) (list)
    """

    motifs = {
        "name": organism_name,
        "ITS_region": [],
        "leader": [],
        "d1d1": [],
        "sp_d2d3_sp": [],
        "tRNA_ile": [],
        "sp_v2_sp": [],
        "tRNA_ala": [],
        "BoxB": [],
        "BoxA": [],
        "D4": [],
        "V3": [],
    }

    # Sequence to be used for the motif search.
    # Found motifs get removed from the seq before moving on to the next one.
    its_seq = slice_its_region(seq_input)
    if its_seq == -1:
        print(
            Back.RED
            + Fore.WHITE
            + "Found ITS Region length is too short (not accurate). Skipping\n"
            + str(organism_name)
        )
        return None

    # store the whole ITS region sequence in the dict. Remains static
    motifs["ITS_region"].append(its_seq)

    d1d1 = get_d1d1("GACCT", r"AGGTC", its_seq)
    if d1d1 is None:
        d1d1 = get_d1d1("GACCA", r"[AT]GGTC", its_seq)
    if d1d1 is None:
        d1d1 = get_d1d1("GACCG", r"[AC]GGTC", its_seq)
    if d1d1 is None:
        d1d1 = get_d1d1("GACCC", r"[AC]GGTC", its_seq)
    if d1d1 is None:
        d1d1 = get_d1d1("GACAA", r"[AT]TGTC", its_seq)  # pico cyano?
    if d1d1 is None:
        d1d1 = get_d1d1("GACTT", r"[AT]GGTC", its_seq)

    if d1d1 is None or len(d1d1) == 0:
        motifs["leader"] = None
        motifs["d1d1"] = None
    else:
        # get the starting index of the d1d1
        leader_start = its_seq.find(d1d1[0])
        # use the above to define where leader ends (start of d1d1)
        motifs["leader"].append(its_seq[0:leader_start])
        for seq in d1d1:
            # Append the d1d1 results to the d1d1 dict key (an array).
            motifs["d1d1"].append(seq)
            # Trim the its_seq. rindex picks the last index of the found string.
            # Equivalent to the end of the d1d1 string.
        its_seq = its_seq[its_seq.rindex(motifs["d1d1"][-1]) :]

    # spacer-D2-spacer D3-spacer region and  tRNA-1
    trna_ile = find_motif("GGGC[TC]ATTA", "GGCCCA", its_seq, 70, 80)
    if trna_ile is None:
        motifs["tRNA_ile"] = None
        motifs["sp_d2d3_sp"] = None
    else:
        # Store what is between the end of d1d1 (first char) and the start of tRNA_ile
        motifs["sp_d2d3_sp"].append(its_seq[0 : its_seq.find(trna_ile[0])])
        for seq in trna_ile:
            motifs["tRNA_ile"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA_ile"][-1]) :]

    # find tRNA-Ala(2)
    trna_ala = find_motif("GGGG", "[TC]CTCCA", its_seq, 70, 80)
    if trna_ala is None:
        motifs["tRNA_ala"] = None
        motifs["sp_v2_sp"] = None
    else:
        # everything between end of tRNA_ile and start of tRNA_ala
        motifs["sp_v2_sp"].append(its_seq[0 : its_seq.find(trna_ala[0])])
        for seq in trna_ala:
            motifs["tRNA_ala"].append(seq)
        its_seq = its_seq[its_seq.rindex(motifs["tRNA_ala"][-1]) :]

    # find BoxB
    boxb = find_motif("CAGC", "GCTG", its_seq, 18, 50)
    if boxb is None:
        boxb = find_motif("AGCA", "CTG", its_seq, 18, 50)
    if boxb is None:
        motifs["BoxB"] = None
    else:
        for seq in boxb:
            motifs["BoxB"].append(seq)
        # rindex picks the last index of the found string.
        # The first BoxB [0] is the longest, so using that.
        its_seq = its_seq[its_seq.rindex(motifs["BoxB"][0]) :]

    # find Box-A
    # reverse the string, look from the back to the front
    boxa = find_motif("TCAAA[GA]G", "A[AC]G", its_seq[::-1], 1, 10)
    if boxa is None:
        motifs["BoxA"] = None
    else:
        for seq in boxa:
            seq = seq[::-1][:-3] #reverse it and trim the ACT that belongs to D4
            motifs["BoxA"].append(seq)  # reverses it back.
        its_seq = its_seq[its_seq.rindex(motifs["BoxA"][0]) - 3 :]

    # find D4
    # TODO; make sure it starts at ACT
    d4 = find_motif("ACT", "TA[TGAC]", its_seq, 4, 20)
    if d4 is None:
        motifs["D4"] = None
    else:
        # D4 is not very accurate, store only the first one.
        motifs["D4"].append(d4[0])
        its_seq = its_seq[its_seq.rindex(motifs["D4"][0]) :]

    # Find v3
    v3 = find_motif("GT[CA]G", "CA[CG]A[GC]", its_seq, 4)
    if v3 is None:
        motifs["V3"] = None
    else:
        # V3 not very accurate, stores only the first one.
        motifs["V3"].append(v3[0])
        its_seq = its_seq[its_seq.rindex(motifs["V3"][0]) :]

    return motifs


def find_motif(start_pattern, end_pattern, seq, min_length=4, max_length=200):

    """Helper function for slice_motifs to find a single motif given basal clamps in the arguments.
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
        start_match = search(start_pattern, seq)
        # If the start position is not found, return None to the main program
        if start_match is None:
            break
        seq = seq[start_match.start() :]
        # Find all the matching bases to the pattern in the argument passed to the function (end)
        end_matches = finditer(r"" + end_pattern, seq)
        # For each match, add the end position to the d1d1_results array as a (start,end) tuple.
        for match in end_matches:
            motif_seq = str(seq[: match.end()])
            if len(motif_seq) > min_length and len(motif_seq) < max_length:
                motif_results.append(motif_seq)
        seq = seq[3:]
    if len(motif_results) == 0:
        return None
    return motif_results


def slice_its_region(seq_input):
    """Helper function that slices the ITS region from the raw sequence extracted from the fasta or genbank data.

    Args:
        seq_input (string): full sequence
        min_length (minimum length): establishes the minimum length that the ITS region can be to avoid false positives

    Returns:
        string: ITS region sequence
    """
    min_length = 20  # Used to filter out bad results.
    its_seq_search = search(r"CCTCCTT", seq_input)
    if its_seq_search is None:
        its_seq_search = search(r"CCTCCTA", seq_input)
    if its_seq_search is None:
        print(
            Fore.YELLOW
            + "\nWARN: Could not find the end of 16S to determine the ITS region boundaries. Results may be inaccurate."
        )
        return seq_input

    if len(seq_input[its_seq_search.start() :]) < min_length:
        print(
            Back.RED + Fore.WHITE + "Region length too short. Skipping this sequence."
        )
        return -1
    return seq_input[its_seq_search.start() + 7 : its_seq_search.start() + 700]


def get_d1d1(start_pattern, end_pattern, sequence):
    """Helper Function to search for d1d1 motif which has different constraints than others.

    Args:
        start (string): d1d1 start pattern
        end (string): d1d1 end pattern
        seq (string): sequence where the d1d1 is to be found.

    Returns:
        list: list containing all the possible d1d1 sequences.
    """
    d1d1_results = []
    # Limits the area in which the d1d1 region is usually found.
    d1d1_search_area = sequence[0:150]
    d1d1_start_position = d1d1_search_area.find(start_pattern, 0, 20)
    if d1d1_start_position == -1:
        return None
    end_matches = finditer(end_pattern, d1d1_search_area)
    # For each match, add the end position to the d1d1_results array as a (start,end) tuple.
    for match in end_matches:
        d1d1_results.append(sequence[d1d1_start_position : match.end()])
    if len(d1d1_results) == 0:
        return None
    return d1d1_results


def print_motifs(motif_list):
    """Gets the dict of motifs and prints them out in the terminal

    Args:
        motifs (dict): list of motifs to print
    """
    if motif_list is None:
        print(Back.RED + Fore.WHITE + "ITS Region not found in this sequence.")
    else:
        print(Fore.WHITE + Style.BRIGHT + motif_list["name"])
        for key in motif_list:
            if key != "name":
                if motif_list[key] is None:  # If the key is empty
                    if key in ("tRNA_ile", "tRNA_ala"):
                        print(
                            Fore.LIGHTCYAN_EX
                            + Style.BRIGHT
                            + key
                            + Fore.RED
                            + " Not present in this operon."
                        )
                        print("\n")
                    else:
                        print(
                            Fore.LIGHTCYAN_EX
                            + Style.BRIGHT
                            + key
                            + Fore.RED
                            + " Not found in this sequence."
                        )
                        print("\n")
                else:
                    if (
                        len(motif_list[key]) > 1
                    ):  # Check if there is more than one sequence in that key
                        print(
                            Fore.CYAN
                            + Style.BRIGHT
                            + key
                            + " \n\t"
                            + Fore.RED
                            + Style.BRIGHT
                            + str(len(motif_list[key]))
                            + " possible sequences found! \n"
                        )
                        for index, entry in enumerate(motif_list[key]):
                            print(
                                Fore.MAGENTA
                                + Style.BRIGHT
                                + "\tSequence "
                                + str(index + 1)
                                + ": "
                                + Fore.LIGHTYELLOW_EX
                                + Style.NORMAL
                                + str(entry)
                                + Style.BRIGHT
                                + Fore.LIGHTMAGENTA_EX
                                + "\n\tLength: "
                                + Style.NORMAL
                                + Fore.LIGHTYELLOW_EX
                                + str(len(entry))
                                + "\n"
                            )
                    else:
                        if key is list(motif_list.keys())[0]:
                            print(Fore.CYAN + Style.BRIGHT + "ITS Region:")
                        else:
                            print(Fore.CYAN + Style.BRIGHT + key + ":")
                        for entry in motif_list[key]:
                            print(
                                Fore.MAGENTA
                                + Style.BRIGHT
                                + "\tSequence "
                                + Fore.LIGHTYELLOW_EX
                                + Style.NORMAL
                                + str(entry)
                                + Style.BRIGHT
                                + Fore.LIGHTMAGENTA_EX
                                + "\n\tLength: "
                                + Style.NORMAL
                                + Fore.LIGHTYELLOW_EX
                                + str(len(entry))
                                + "\n"
                            )
    terminal_width = get_terminal_size()[0]
    print(Fore.LIGHTGREEN_EX)
    for x in range(terminal_width):
        print(Fore.LIGHTGREEN_EX + "-", end="")
    print("\n\n")


def trna_check(motifs):
    """Checks for presence of tRNAs to check for homologous operons.

    Args:
        motifs (dict): list of motifs of one organism

    Returns:
        dict: dictionary with the tRNAs and if they are found or not.
    """
    trnas = {"tRNA_ile": "Not Found", "tRNA_ala": "Not Found"}
    if motifs["tRNA_ile"] is not None:
        trnas["tRNA_ile"] = "Found"
    if motifs["tRNA_ala"] is not None:
        trnas["tRNA_ala"] = "Found"
    print(Fore.CYAN + "\n" + str(list(motifs.keys())[0]) + "\n")
    if trnas["tRNA_ile"] == "Found":
        print(Fore.GREEN + "tRNA1: " + trnas["tRNA_ile"])
    else:
        print(Fore.RED + "tRNA1: " + trnas["tRNA_ile"])
    if trnas["tRNA_ala"] == "Found":
        print(Fore.GREEN + "tRNA1: " + trnas["tRNA_ile"])
    else:
        print(Fore.RED + "tRNA1: " + trnas["tRNA_ile"])


def generate_fasta(organisms, bymotif=False):
    """Writes output to fasta file

    Args:
        organisms (list): list of organisms to write
        singlefile (bool, optional): If -s is enabled, output to one file per motif with multiple organisms. Defaults to False.
    """

    p = Path("output")
    # If ./output doesn't exist, create it.
    p.mkdir(exist_ok=True)
    for org in organisms:
        if not bymotif:
            fasta_output = open(
                "./output/" + org["name"].replace(".", "") + "-motifs.fasta",
                "w",
                encoding="utf8",
            )
        for key in org:
            if key != "name":
                if bymotif:
                    fasta_output = open(
                        "./output/" + key + "-output.fasta", "a+", encoding="utf8"
                    )
                    fasta_output.write('\n')
                if org[key] is not None:
                    if len(org[key]) > 1:
                        for i in range(len(org[key])):
                            fasta_output.write(
                                ">"
                                + org["name"]
                                + " - "
                                + str(key)
                                + "."
                                + str(i)
                                + "\n"
                                + str(org[key][i])
                                + "\n"
                            )
                    else:
                        for entry in org[key]:
                            fasta_output.write(
                                ">"
                                + org["name"]
                                + " - "
                                + str(key)
                                + "\n"
                                + str(entry)
                                + "\n"
                            )

        fasta_output.close()


# Main
# Create the parser
parser = argparse.ArgumentParser(
    description="Find motifs within Cyanobacteria ITS region sequences"
)
# Create a group of mutually exclusive argument options
group = parser.add_mutually_exclusive_group()
group.add_argument(
    "-f",
    "--fasta",
    action="store",
    help="Find motifs in a given fasta file.",
)
group.add_argument(
    "-g",
    "--genbank",
    help="Fetch a sequence from a given genbank accession number.",
    nargs="+",
)
parser.add_argument(
    "-s",
    "--select",
    help="Select which motifs to extract",
    default="all",
    nargs="*",
    choices=(
        "ITS_region",
        "leader",
        "d1d1",
        "sp_d2d3_sp",
        "tRNA_ile",
        "sp_v2_sp",
        "tRNA_ala",
        "BoxB",
        "BoxA",
        "D4",
        "V3",
        "all",
    ),
)
parser.add_argument(
    "-e",
    "--email",
    help="Provide email for genbank query. If not provided, you will be prompted for one at runtime.",
    default=None,
)
parser.add_argument(
    "-t",
    "--trna",
    help="Returns ONLY how many tRNAs are found to use for homologous operon verification.",
    action="store_true",
)
parser.add_argument(
    "-o", "--output", help="Outputs a Fasta file with the results", action="store_true"
)

# Process the arguments
# If there is no argument then print the help
if len(argv) == 1:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

if args.fasta:
    try:
        fasta = open(args.fasta, "r", encoding="UTF-8")
        sequences = parse_fasta(fasta)
        processed_organisms = process_seqs(sequences)
    except IOError:
        system_exit("File not found. Exiting")

if args.genbank:
    # Get email for Entrez. Either from argument or ask the user for input
    if args.email:
        user_email = args.email
    else:
        print(
            "To query the Entrez database via a script, NCBI requires attaching an email address to the query."
        )
        print(
            "This script does not check if the email is valid, does not store it and does not share it with anyone else."
        )
        print(
            "The Biopython library is used to handle the communication with Entrez and the handling of the provided email address"
        )
        user_email = input("Please enter your email for the Entrez query: ".lower())
    while not valid_email(user_email):
        print("Invalid email, try again: ")
        user_email = input(
            "Entrez requires an email to be associated with the requests. Please enter your email: ".lower()
        )
    sequences = []
    for gb in args.genbank:
        # Required to not go over the 3 queries/second threshold imposed by Entrez
        sleep(0.35)
        try:
            sequences.append(parse_genbank(gb, user_email))
            processed_organisms = process_seqs(sequences)
        except IOError:
            system_exit("Error while parsing accession number. Exiting.")

if args.trna:
    for organism in processed_organisms:
        trna_check(organism)
else:
    if args.select:
        delete_items = [
            "ITS_region",
            "leader",
            "d1d1",
            "sp_d2d3_sp",
            "tRNA_ile",
            "sp_v2_sp",
            "tRNA_ala",
            "BoxB",
            "BoxA",
            "D4",
            "V3",
        ]
        if args.select == "all":
            delete_items.clear()
        else:
            delete_items = [x for x in delete_items if x not in args.select]
        for organism in processed_organisms:
            for item in delete_items:
                del organism[item]
            print_motifs(organism)
    if args.output:
        generate_fasta(processed_organisms, args.select != "all")
