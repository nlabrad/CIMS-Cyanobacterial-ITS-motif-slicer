# The Motif Deli (?)

Pre-Requirements:

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

If running python --version returns Python 3.X.X:

pip install colorama
pip install biopython

## Windows:

[Follow this to install pip](https://www.liquidweb.com/kb/install-pip-windows/)

### Install the dependencies:

pip install colorama
pip install biopython

Navigate to the location of the file:
For example:

cd C:/Users/{your-username}/Desktop



## Running the app:
(For MacOS, you may have to run these using python3 instead of python)

For help:

python motifsearch.py -h or python motifsearch.py --help
### Processing a fasta file:

python motifsearch.py -f {filename}

For example:

python motifsearch.py -f microbes.fasta

### Fetching a file from Genbank then processing

python motifsearch.py -g {accession Number}

For example:

python motifsearch.py -g MT425922.1


If you run into any issues please open an issue ticket here: https://github.com/nlabrad/motifsearch/issues

python motifsearch.py -g


