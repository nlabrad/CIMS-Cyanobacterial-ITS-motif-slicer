
<div align='center'>
  
# CIMS: Cyanobacterial ITS Motif Slicer

CIMS is a tool to extract the commonly used ITS folding motifs from a 16s-23s rRNA sequence. It takes a fasta or at least one Genbank accession number and returns a list of motifs with their respective labels for each of the sequences provided. 
Dedicated to the cyanobacteria researches who spend many hours highlighting motifs in MS Word.

</div>

### Table of content
========
 * [Why did we make CIMS?](#why)
 * [What does it do again?](#what)
 * [Installation](#installation)
 * [Usage](#usage)
 * [Possible Errors](#possible-errors)



### Why did we make this tool?

The 16S-23S rRNA internal transcribed spacer (ITS) is a commonly employed phylogenetic marker in cyanobacterial systematics. Examination of ITS regions allows researchers to discover congruencies and apomorphies between species of cyanobacteria. This gives the researcher more evidence when erecting new cryptic taxon or analyzing previously unresolved taxonomic relationships. The challenge however is that historically researchers must manually dig through sequence data to visually find and identify ITS sequence motifs. This painstaking process deters researchers from using ITS motifs, leads to errors, and not to mention… causes headaches.

We knew there was a better way to do this, so after dissecting the manual process, we created *CIMS*. 
CIMS finds the commonly used ITS folding motifs such as D1-D1’, Box B, tRNA-ile and tRNA-ala to ensure researchers are using homologous operons when comparing ITS secondary structures between taxa. 


### What does it do again?

+ CIMS is a terminal application written in Python that 
+ CIMS can process one or more Genbank accession numbers or a fasta file with one or more sequences. 
+ It automatically talks to Genbank for you so you don't have to download the fasta files yourself.
+ Returns a text output with the motifs identified and their lenghts for you to use as you please.

In the current version of the software, the motifs included in the standard output are:
  +	Leader
  +	D1-D1`
  +	Spacer – D2 – Spacer 
  +	tRNA-ala
  +	Spacer – V2 - Spacer
  +	tRNA- ile 
  +	Box B 
  +	D4
  +	BoxA 
  +	V3 
                  
### Installation


  #### Pre-Requisites
    We get it, you're a biologist, we got you. All you need is beginner level of terminal... maybe not even that much. If you know how to browse to a directory (```cd```) and run an executable(```./cims```), you're good to go.


  #### Simple Method: Download the pre-packaged files from [Releases](https://github.com/nlabrad/CIMS-Cyanobacterial-ITS-motif-slicer/releases).

    To keep things simple, we pre-packaged ```CIMS``` with all it's dependencies into a single file and compiled it for Windows, Linux and MacOS. These files are available under [Releases](https://github.com/nlabrad/CIMS-Cyanobacterial-ITS-motif-slicer/releases).

    + Download the zip file that corresponds to your system.
    + Unzip it in whichever directory you'd like. 
    + You're done! :open_hands:
    + To run, open your favorite Terminal, ```cd``` to that directory, and run cims as an executable, usually by typing ```./cims```.

    To keep things simple, we suggest saving cims to the directory where you'll have the fasta files you want to process. 
    If you're pulling your sequences straight from Genbank, it doesn't really matter.

  #### Advanced Method: Download the Python script.

    If you want to perhaps make your own changes to the flanking regions, or make changes to the code, you can simply download CIMS.py from and run it with Python. (But you probably already knew that if that's what you wanted). 

    To run ```CIMS``` you will need:
    + Python 3
    + BioPython: ```$ pip install Biopython```
    + Colorama ```$ pip install colorama```

    [BioPython](https://biopython.org/) allows CIMS to communicate with Genbank to download sequences.
    [Colorama](https://github.com/tartley/colorama) allows us to easily output the motifs in pretty colors.

  Once you have those dependencies installed (either globally or in a virtual environment), simply run cims.py.

### Usage

  CIMS runs in the terminal. It is provided a sequence either through a FASTA file or by fetching them from Genbank based on accession numbers.
  The input for this tool must either be a fasta file with one or more properly formatted 16s-23s ITS sequences or a Genbank accession number to a 16s-23s ITS sequence.


  Navigate to the location of the file (either cims or cims.py if you downloaded the script):

  For example, in Windows, you'd use ```cd``` to move to a directory as such:

  ```cd C:/Users/{your-username}/Desktop/PathtoFile```

  Or in Linux/Mac:

  ```cd /home/{your username}/{where you downloaded cims}```

  To run CIMS, simply execute it by running ```./cims``` or ```python cims.py``` from the directory where it was saved. 

  When running this on your terminal the output will include all motifs found in the sequences given to the program. If you would like to save the output of your run remember to use “>>” to save output into a text file:

  ``` cims -f myfasta.fasta >> motifs.txt``` 

  The list of flags, arguments and their descriptions are below:

  ```shell
  Usage: cims [OPTIONS]

  Options:
  -f, --fasta PATH-TO-FASTA-FILE                                             Provide FASTA to be processed.
  -g, --genbank ACCESSION1 [ACCESSION2 ...]                                  Provide one or more Genbank Accession Numbers to fetch and process.
  -s, --select {{leader,d1d1,sp_v2_sp,trna_ile,trna_ala,boxa,boxb,d4,v3,all} Select which motifs to print out. By default it prints all.
  -e, --email                                                                Provide an email to be used when querying Genbank. An NCBI requirement.
  -j, --json                                                                 Create a json file in the working directory with the output.
  -t, --trna                                                                 Returns ONLY how many tRNAs were found per sequence. 
  ```

  #### Examples:
  ```cims =f allmycyanos.fasta```

  Result: CIMS will process the provided fasta file and return all the motifs it finds.

  ```cims -f ~/home/me/fasta/limnothrix_16-23_ITS.fasta -s d1d1, trna_ile, trna_ala, boxb```

  Result: Processes the limnothrix_16-23_ITS.fasta file stored in a directory that resides in /home/me/fasta and asks CIMS to only output d1d1, the tRNAs and BoxB motifs.

  ```cims -g KU574618.1 -e my@email.com```

  Result: Fetches the sequence of KU574618.1 from Genbank (providing an email that is required by NCBI), processes the sequence, and returns the motifs.

  ```cims -f allmycyanos.fasta -t```

  Result: Fetches the sequence from Genbank, and returns how many tRNAs were found on each organism. This allows to easily check if the organisms in the fasta are homologous operons.

  >If you ever get lost, you can always run ```cims -h``` or ```python cims.py -h``` and you will get a quick reference of the available options.

### Possible errors: 

  #### 1. ```“Could not find the end of 16S to determine the ITS region boundaries”```
  This error means that the sequence given to the software did not contain the sequence that represents the end of the 16S region (CCTCCTT). You may proceed with the run if you have fed the program the ITS region only and everything will run as normal otherwise, abort the run for that sequence by typing “N” when prompted “Proceed with search anyway? (Y/N)”. This will allow the program to move onto the next sequence in the fasta file or allow you to try again with another file/accession #. 

  #### 3. ```“Region length too short. Skipped.”```
  This will be printed if the ITS region after the end of the 16S gene is under 20bps. This feature is coded to remove sequences with ITS regions that are too small to be used to find any of the motifs. 

  #### 4. ```“Not found in this sequence.” ```
  This output will be printed when a particular motif was not found in the ITS sequence. This could be because the flanking regions are unique or otherwise rare and so the software did not find these. If this happens frequently in your dataset, please report this to us in the “Issues” page of the GitHub so that we can address this error and improve the code.

  #### 5. ```“Not present in this operon” ```
  This will be printed only regarding tRNAs in the sequence. If the program does not find tRNA-ala or tRNA-ile, it will assume that this operon        does not contain one or both tRNAs. Remember, it is best to use homologous operons when comparing ITS motifs between taxa (ie. Operons containing        the same number of tRNAs). 


