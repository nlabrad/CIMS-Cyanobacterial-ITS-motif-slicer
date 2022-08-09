import re

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
    #returns an array of tuples. TODO: any better options?
    
def search_motif(seq_start, seq_end, pattern, seq):
    
    search_results = []
    # Limits the area in which the d1d1 region is found, usually.
    search_area = seq[seq_start:seq_end]


    start_matches = re.finditer(pattern, search_area)
    for match in start_matches:
        search_results.append(seq[match.start():match.end()])
    #If the start position is not found, return None to the main program
    if (start_matches == None):
        return None
    # But if it is found, try searching for the end pattern.
    
    return search_results
    #returns an array of tuples. TODO: any better options?