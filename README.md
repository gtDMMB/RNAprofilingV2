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

This code has required dependencies on the numpy, networkx, matplotlib, and pygraphviz
python libraries. It has been tested with numpy 1.19.5, networkx 2.4, matplotlib 3.3.4, 
pygraphviz 1.6, and python 3.6.9. It also requires an installation of graphviz
to generate the output diagrams.

There is an optional dependency on the RNAlib (Vienna)
python bindings. If the are not installed, the code requires either a pregenerated 
sample or an install of RNAstructure in order to run. RNAlib 2.5.1 or newer 
must be installed for seeds to work while using its sampling. 

To run this code without installing it locally visit https://rnaprofiling.gatech.edu/.

Citation information coming soon.

