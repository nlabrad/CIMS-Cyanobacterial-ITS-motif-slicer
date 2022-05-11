import re
from Bio import Entrez, __version__, SeqIO, Seq #import the required libraries


def fetch(accession):
    Entrez.email = 'nlab@fastmail.com' # email reported to entrez to associate with the query
    with Entrez.efetch(db='nucleotide', id=accession, rettype='gb') as handle: # id is what genbank ID to query, type is genbank
        record = SeqIO.read(handle, "gb") 
        genbank_organism_taxonomy = record.annotations["taxonomy"]
        #genbank_organism_sequence = record.seq
        print ("Taxonomy: ", genbank_organism_taxonomy)
        # print ("Sequence:", genbank_organism_sequence)

def get_tax_id(organism):
    Entrez.email = 'nlab@fastmail.com' # email reported to entrez to associate with the query
    organism = organism.replace(" ", "+").strip() # replace whitespaces with +
    handle = Entrez.esearch(term = organism, db = "taxonomy", retmode = "xml") #search for the organism
    record = Entrez.read(handle)
    return record['IdList'][0] #IdList[0] contains the taxID

def get_taxonomy(tax_id):
    with Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml') as handle:
        record = Entrez.read(handle, validate=False)
        print(record)
        return record
        

#input: fasta file or search genbank
#output: motifs
# 1 Find the ITS region first. 16S ends with CCTCCTT. Ends in CAGAC: 16S...CCTCCTT-ITS...CAGAC-23S...
# 2 Find motifs in ITS region

#sample string with the 16S-23S region
seq = "AGTTGATCCTGGCTCAGGATGAACGCTGGCGGTATGCTTAACACATGCAAGTCGAACGGTCTCTTCGGAGATAGTGGCGGACGGGTGAGTAACGCGTGAGAATCTAGCTTCAGGTTCGGGACAACCACTGGAAACGGTGGCTAATACCGAATGTGCCGAGAGGTGAAAGGCTTGCTGCCTGAAGATGAGCTCGCGTCTGATTAGCTAGTTGGTGGGGTAAAAGCCTACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATACCGCGTGAGGGAGGAAGGCTCTTGGGTCGTAAACCTCTTTTCTCAGGGAAGAAAAAAATGACGGTACCTGAGGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCCGCAGGTGGCGATGTAAGTCTGCTGTTAAAGAGTGAGGCTCAACCTCATAAAAGCAGTGGAAACTACATGGCTAGAGTGCGTTCGGGGTAGAGGGAATTCCTGGTGTAGCGGTGAAATGCGTAGAGATCAGGAAGAACACCGGTGGCGAAAGCGCTCTGCTAGGCCGCAACTGACACTGAGGGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTAGCCGTAAACGATGGATACTAGGTGTGGCTTGTATCGACCCGAGCCGTGCCGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCAAGACTTGACATGTCGCGAATCCCGGTGAAAGCTGGGAGTGCCTTCGGGAGCGCGAACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTTTTTAGTTGCCAGCATTAAGTTGGGCACTCTAGAGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCAGCATGCCCCTTACGTCTTGGGCTACACACGTACTACAATGCTACGGACAAAGGGCAGCTACACAGCGATGTGATGCAAATCCAAGAAACCGTAGCTCAGTTCAGATCGCAGGCTGCAACTCGCCTGCGTGAAGGAGGAATCGCTAGTAATTGCAGGTCAGCATACTGCAGTGAATTCGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGAAGCTGGTCGCGCCCGAAGTCATTACCCCAACCTTTCGAGGAGGGGGATGCCTAAGGCAGGACTGGTGACTGGGGTGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGTGGCTGGATCACCTCCTTTTTAGGGAGACCTACACCCCTCAAATCTTGAAAGCAAAATGCAAATAGAGATTGAGTTGGTCTATCCTAGGTCGGTCGCAGATAATTGTTGAAGCTTTCAAACTATGATTTGGTTCGATATGGGCTATTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCCCTGGTTCGAGTCCAGGATGGCCCACCTGAAAGTAGTTAGTAATTGCTAATTTATAATTCGTAAATAACTACGAATTAAGGATTAAGAATTACGAATTATAATCTGGGGGTTTAGCTCAGTTGGTAGAGCGCCTGCTTTGCAAGCAGGATGTCAGCGGTTCGAGTCCGCTAACCTCCACATTGGAAAAAAGCCAGCAAAGTGGTGAAAACCAGACTGCTGGGAGAATAATCCAGCCAGAACCTTGAAAACTGCATAGTAACGCGAAAAATAGCAGGCAGACACAGACATCGAAAAGATGAAAGTGAAAGCAGGTGAAAACCAATGAAATTGTGGTCAAGCTAATAAGGGCTAACGGTGGATACCTAGGCACACAGAG"
#seq = "AGTTGATCCTGGCTCAGGATGAACGCTGGCGGTATGCTTAACACATGCAAGTCGAACGGTCTCTTCGGAGATAGTGGCGGACGGGTGAGTAACGCGTGAGAATCTAGCTTCAGGTTCGGGACAACCACTGGAAACGGTGGCTAATACCGAATGTGCCGAGAGGTGAAAGGCTTGCTGCCTGAAGATGAGCTCGCGTCTGATTAGCTAGTTGGTGGGGTAAAAGCCTACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATACCGCGTGAGGGAGGAAGGCTCTTGGGTCGTAAACCTCTTTTCTCAGGGAAGAAAAAAATGACGGTACCTGAGGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCCGCAGGTGGCGATGTAAGTCTGCTGTTAAAGAGTGAGGCTCAACCTCATAAAAGCAGTGGAAACTACATGGCTAGAGTGCGTTCGGGGTAGAGGGAATTCCTGGTGTAGCGGTGAAATGCGTAGAGATCAGGAAGAACACCGGTGGCGAAAGCGCTCTGCTAGGCCGCAACTGACACTGAGGGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTAGCCGTAAACGATGGATACTAGGTGTGGCTTGTATCGACCCGAGCCGTGCCGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCAAGACTTGACATGTCGCGAATCCCGGTGAAAGCTGGGAGTGCCTTCGGGAGCGCGAACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTTTTTAGTTGCCAGCATTAAGTTGGGCACTCTAGAGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCAGCATGCCCCTTACGTCTTGGGCTACACACGTACTACAATGCTACGGACAAAGGGCAGCTACACAGCGATGTGATGCAAATCCAAGAAACCGTAGCTCAGTTCAGATCGCAGGCTGCAACTCGCCTGCGTGAAGGAGGAATCGCTAGTAATTGCAGGTCAGCATACTGCAGTGAATTCGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGAAGCTGGTCGCGCCCGAAGTCATTACCCCAACCTTTCGAGGAGGGGGATGCCTAAGGCAGGACTGGTGACTGGGGTGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGTGGCTGGATCACCTCCTTTTTAGGGAGACCCACACCCCTCAAATCTTGAAAGCAAAATGCAAATAGAGATTGAGTTGGTCTATCCTAGGTCGGTCGCAGATAATTGTTGAAGCTTTCAAACTATGATTTGGTTCGATATGGGCTATTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCCCTGGTTCGAGTCCAGGATGGCCCACCTGAAAGTAGTTAGTAATTGCTAATTTATAATTCGTAAATAACTACGAATTAAGGATTAAGAATTACGAATTATAATCTGGGGGTTTAGCTCAGTTGGTAGAGCGCCTGCTTTGCAAGCAGGATGTCAGCGGTTCGAGTCCGCTAACCTCCACATTGGAAAAAAGCCAGCAAAGTGGTGAAAACCAGACTGCTGGGAGAATAATCCAGCCAGAACCTTGAAAACTGCATAGTAACGCGAAAAATAGCAGGCAGACACAGACATCGAAAAGATGAAAGTGAAAGCAGGTGAAAACCAATGAAATTGTGGTCAAGCTAATAAGGGCTAACGGTGGATACCTAGGCACACAGAG"
#Extraction of ITS region from a 16S-23S region.

# Find CCTCCTT, get rid of what's before that.
def extractITS(geneRegion):
    if(re.search("CCTCCTT", geneRegion) == None): # If It cannot find the start of the ITS region4
        print("ITS Region not found in this sequence")
        return None
    else:    
        itsStartPosition = re.search("CCTCCTT", geneRegion).start() + 7
        itsRegion = geneRegion[itsStartPosition:] #Store everything from CCTCCTT to the end
        #itsEndSearch = re.search("CAGAC", itsRegion) #Search for CAGAC, the end of the ITS
        #if (itsEndSearch != None):
            #itsRegion = itsRegion[0:itsEndSearch.start() + 5] #Region = what found before including the CAGAC (+5 from start)
        #If itsEndSearch is none, then include everything after CCTCCTT in the region.
        print("\nThe ITS region is:", itsRegion)
        print("ITS region length", len(itsRegion))
        return itsRegion

itsregion = extractITS(seq)

def extractMotifs(its_region):
    
    motifs = {} #dictionary. motifs[key]=[sequence,length]
         
    
    d1d1Search = re.search(r"GACC[CT](.*?)AGGTC",its_region) #find text between basal clamps, starting with GACCT/C to the first AGGTC (*? is lazy search)
    if (d1d1Search == None):
        d1d1Search = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (d1d1Search == None):
        motifs["d1d1"] = None
    else: 
        motifs["leader"] = [its_region[0:d1d1Search.start()], len(its_region[0:d1d1Search.start():])]
        motifs["d1d1"] = [its_region[d1d1Search.start():d1d1Search.end()],len(its_region[d1d1Search.start():d1d1Search.end()])]
        its_region = its_region[d1d1Search.end():] #trim processed its region
        print("\nD1D1: ",motifs["d1d1"])
        print("D1D1 lenght: ", motifs["d1d1"][1])
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))
    
    # spacer-D2-spacer D3-spacer region
    
    # find tRNA-Ile(1) and tRNA-Ala(2)
    tRNA1Search = re.search(r"GGGCTATTA(.*?)GGCCCA",its_region) #find text between basal clamps
    if (tRNA1Search == None):
        tRNA1Search = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (tRNA1Search == None):
        motifs["tRNA1"] = None
    else:  
        motifs["d2d3"] = [its_region[0:tRNA1Search.start()], len(its_region[0:tRNA1Search.start():])]
        print("\nspacer-D2-spacer D3-spacer", motifs["d2d3"])
        print("spacer-D2-spacer D3-spacer lenght:", motifs["d2d3"][1])
        motifs["tRNA1"] = [its_region[tRNA1Search.start():tRNA1Search.end()],len(its_region[tRNA1Search.start():tRNA1Search.end()])]
        print("\ntRNA1: ",motifs["tRNA1"])
        print("tRNA1 lenght: ", motifs["tRNA1"][1])
        
        its_region = its_region[tRNA1Search.end():] #trim processed its region
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))
    
    
    tRNA2Search = re.search(r"GGGG(.*?)CCTCCA",its_region)  #find text between basal clamps
    if (tRNA2Search == None):
        tRNA2Search = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (tRNA2Search == None):
        motifs["tRNA2"] = None
    else:  
        motifs["tRNA2"] = [its_region[tRNA2Search.start():tRNA2Search.end()], len(its_region[tRNA2Search.start():tRNA2Search.end()])] #motif[] = [seq,len]
        print("\ntRNA2: ",motifs["tRNA2"])
        print("tRNA2 lenght: ", motifs["tRNA2"][1])
        its_region = its_region[tRNA2Search.end():] #trim processed its region
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))
        
    #find BoxB 
    BoxBSearch = re.search(r"AGCA(.*?)TGCT", its_region) #find text between basal clamps
    if (BoxBSearch == None):
        BoxBSearch = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (BoxBSearch == None):
        motifs["BoxB"] = None
    else:  
        motifs["BoxB"] = [its_region[BoxBSearch.start():BoxBSearch.end()], len(its_region[BoxBSearch.start():BoxBSearch.end()])] #motif[] = [seq,len]
        print("\nBoxB: ",motifs["BoxB"])
        print("BoxB lenght: ", motifs["BoxB"][1])
        its_region = its_region[BoxBSearch.end():] #trim remaining its region
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))    
        
    #find Box-A and D4
    BoxASearch = re.search(r"G[CA]A(.*?)GAAA",its_region) #find text between basal clamps
    if (BoxASearch == None):
        BoxASearch = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (BoxASearch == None):
        motifs["BoxA"] = None
    else:  
        motifs["BoxA"] = [its_region[BoxASearch.start():BoxASearch.end()], len(its_region[BoxASearch.start():BoxASearch.end()])] #motif[] = [seq,len]
        print("\nBoxA: ",motifs["BoxA"])
        print("BoxA lenght: ", motifs["BoxA"][1])
        its_region = its_region[BoxASearch.end():] #trim remaining its region
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))    
        
        
    D4Search = re.search(r"ACT(.*?)AAA",its_region)
    if (D4Search == None):
        D4Search = re.search(r"GACCA(.*?)[TA]GGTC",its_region) #find text starting with GACCA to the first T/AGGTC (*? is lazy search)
    if (D4Search == None):
        motifs["D4"] = None
    else:  
        motifs["D4"] = [its_region[D4Search.start():D4Search.end()], len(its_region[D4Search.start():D4Search.end()])] #motif[] = [seq,len]
        print("\nD4: ",motifs["D4"])
        print("D4 lenght: ", motifs["D4"][1])
        its_region = its_region[D4Search.end():] #trim remaining its region
        print("\nRemaining ITS Region to be processed:", its_region)
        print("Remaining Length", len(its_region))    
    # end of extractMotif
    
extractMotifs(itsregion)





