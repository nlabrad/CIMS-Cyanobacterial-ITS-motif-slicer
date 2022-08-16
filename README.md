# CIMS: Cyanobacterial ITS motif slicer 

## What are ITS motifs used for? 

The 16S-23S rRNA internal transcribed spacer (ITS) is a commonly employed phylogenetic marker in cyanobacterial systematics. Examination of ITS regions allows researchers to discover congruencies and apomorphies between species of cyanobacteria. This gives the researcher more evidence when erecting new cryptic taxon or analyzing previously unresolved taxonomic relationships. The challenge however is that historically researchers must manually dig through sequence data to visually find and identify ITS sequence motifs. This painstaking process deters researchers from using ITS motifs, leads to errors, and not to mention… causes headaches.


## What does this software do?
This tool has been created to skim through and find the commonly used ITS folding motifs such as D1D1’ and Box B and both 16s-23s ITS tRNAs to ensure researchers are using homologous operons when comparing ITS secondary structures between taxa. 

   ### The program will return the entire ITS sequence and a list of ITS motifs, the corresponding motif sequences, and the motif lengths. In the current version of the software, the motifs included in the standard output are:
                  •	Leader
                  •	D1D1`
                  •	tRNA-ala
                  •	tRNA- ile 
                  •	Box B 
                  
   ### The "full" output will also include the motifs below however as of now finding these motifs is less reliable so use at your own discretion. 
                  •	Spacer – D2 – Spacer 
                  •	Spacer – V2 - Spacer
                  •	D4
                  •	BoxA 
                  •	V3 

## How to use this software: 
The input for this tool must either be a fasta file with one or more properly formatted 16s-23s ITS sequences or a Genbank accession number to a 16s-23s ITS sequence. When running this on your terminal the output will include all motifs found in the sequences given to the program. If you would like to save the output of your run remember to use “>>” to save output into a text file. 


### Running the app:
(For MacOS, you may have to run these using python3 instead of python)

Navigate to the location of the file (motifsearch.py):
For example:

                  cd C:/Users/{your-username}/Desktop/PathtoFile

For help:

                  python3 motifsearch.py -h  
                  python3 motifsearch.py --help

### Processing a fasta file:

                  python3 motifsearch.py -f {filename}

For example:

                  python3 motifsearch.py -f microbes.fasta

### Fetching a file from Genbank then processing

                  python3 motifsearch.py -g {accession Number}

For example:

                  python3 motifsearch.py -g MT425922.1


## Pre-Requirements:

Install Python 3 in whatever environment you're in.

In MacOS, the default Python version is Python 2.7 but it's old and won't work. You still need to download Python 3.

Install in Mac: https://lmgtfy.app/?q=How+to+install+python3+on+mac

Install in Windows: https://lmgtfy.app/?q=How+to+install+python+3+on+windows

## Install dependencies

### MacOS/Linux:

Open a terminal.
If running python --version in your terminal returns Python 2.7 and you have Python 3 installed:

                  pip3 install colorama
                  pip3 install biopython

If running python3 --version returns Python 3.X.X:

                  pip install colorama
                  pip install biopython

## Windows:

[Follow this to install pip](https://www.liquidweb.com/kb/install-pip-windows/)

### Install the dependencies:

                  pip install colorama
                  pip install biopython

If you run into any issues please open an issue ticket here: https://github.com/nlabrad/motifsearch/issues

                  python3 motifsearch.py -g

## Possible errors: 

### 1. “Could not find the end of 16S to determine the ITS region boundaries”  ```
This error means that the sequence given to the software did not contain the sequence that represents the end of the 16S region (CCTCCTT). You        may proceed with the run if you have fed the program the ITS region only and everything will run as normal otherwise, abort the run for that            sequence by typing “N” when prompted “Proceed with search anyway? (Y/N)”. This will allow the program to move onto the next sequence in the fasta        file or allow you to try again with another file/accession #. 

### 3. “Region length too short. Skipped.”
This will be printed if the ITS region after the end of the 16S gene is under 20bps. This feature is coded to remove sequences with ITS              regions that are too small to be used to find any of the motifs. 

### 4. “Not found in this sequence.” 
This output will be printed when a particular motif was not found in the ITS sequence. This could be because the flanking regions are unique or      otherwise rare and so the software did not find these. If this happens frequently in your dataset, please report this to us in the “Issues” page        of the GitHub so that we can address this error and improve the code.

### 5. “Not present in this operon” 
This will be printed only regarding tRNAs in the sequence. If the program does not find tRNA-ala or tRNA-ile, it will assume that this operon        does not contain one or both tRNAs. Remember, it is best to use homologous operons when comparing ITS motifs between taxa (ie. Operons containing        the same number of tRNAs). 
