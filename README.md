# RNAprofiling V2

RNAprofiling takes a sample of RNA secondary structures and 
produces an interactive webpage summarizing them. 

Running tree.py with appropriate options will produce an output\ folder
which contains an index.html file and all necessary auxiliary files. This 
folder can be distributed to anyone to view without any software other 
than an internet browser.

Users may provide a sample in dot or concatenated ct format. If RNAstructure
or ViennaRNA with Python bindings are available locally, the script can 
use them to generate a sample.

If ViennaRNA is available, you may run 

`python3 tree.py sequence.fasta`.

For RNAstructure sampling use 

`python3 tree.py --RNAstructure_location /path/to/rnastructure sequence.fasta`.

To provide a pregenerated sample in dot format 

`python3 tree.py --sample_file /path/to/sample.dot --sample_format dot sequence.fasta`.

To run this code without installing it locally visit https://rnaprofiling.gatech.edu/.

Citation information coming soon.

