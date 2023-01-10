# HelixAnalysis

This repo contains prototype code for generating a decision tree and placing 
it in an interactive webpage. 

Running tree.py produces all of the data needed to populate the website in the 
test\_site folder. It produces an arc diagram folder and a radial diagram 
folder, both of which should be used to replace the folders of the same names 
in the test\_site folder. tree.py also produces 4 javascript files (.js) which 
should be moved to the test\_site/Data folder. There are 3 other files in 
the javascript folder which also need to be updated. At the moment these are 
not automatically generated. One of them contains the name of the sequence, 
one contains the decision tree as svg text and the third contains the entire 
sample in the .gtboltz format. To work around limitations with webpages being 
run locally, all of these are javascript files which have the data stored in 
a variable as a string. 

tree.py also produces a graphviz file with the .dot extension, which should 
be used to generate the svg which goes in the treeSVG.js file in the 
test\_site/Data folder. 

The Decision\_Tree\_Example\_Sites repository contains 5 populated websites 
for various sequences.

