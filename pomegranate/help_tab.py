from dash import dcc

def help_text():
    return dcc.Markdown(''' 

## Need Help?
---
### Welcome to POMEGRANATE!  

POMEGRANATE is an interactive protein structural motif explorer which
allows you to do many things including: 

> - View an abstracted plot of the structure around a phosphorylation site
> - Generate structural motifs of input proteins around a given phosphorylation site
> - View structural motifs in colour and greyscale with interactive information  
> - Compare two or more structural motifs 
> - Explore clustered groupings of structural motifs based on a novel similarity metric (graph neural network  
embeddings) for both labelled (known kinase) and unlabelled (unknown kinase) data

&nbsp

With additional functionality provided through a command-line interface:

> - Retrieve PDB files from AlphaFold V2 given a list of Swiss-Prot IDs 
> - Construct amino-acid graphs based on a thresholded radius from phosphosite, and RSA (residue solvent 
accessibility) cutoff
> - Perform unsupervised representation learning of protein graphs using a Graph Neural Network (GNN)
> - Generate and plot embeddings as image files after performing dimensionality reduction (tSNE, UMAP) 


***

## 1. Page Elements

---

### 1.1 Side Panel

The side panel hosts the input selection for the protein (see 2.1).

---

### 1.2 Different Tabs

The **'Visualisation' tab** is the main screen for the tool where you are able to enter the other elements of the  
required input and view the structural motif generated from the protein (see 2 & 3).

&nbsp 

The **'Compare Motifs' tab** allows you to view 1 or more motifs from the same protein in a side-by-side layout  
so that you can visually compare the similarities and differences (see 4). 

&nbsp 

The **'Proteome View' tab** allows you to view and explore pre-computed clusterings of similar structural  
motifs (see 5). 

&nbsp 

This is the **'Help' tab** which provides documentation and instructions for POMEGRANATE.


*** 

## 2. Input


There are 4 different elements of input required to generate the structural motif: the protein, the specific  
phosphorylation site (phosphosite) on that protein, the radius of interest around that site  
and the Residue Solvent Accessibility threshold. 

---

### 2.1 Protein 

In order to generate a structural motif of a phosphorylation site, you will need to choose the protein that you are  
interested in.  

&nbsp  

In order to select the protein from the database, enter the ID of the protein that you wish to select. The tool will  
then retrieve the file of the protein you specified from the relevant database. 

---

### 2.2 Phosphosite 

There are many options of how to select the phosphosite of interest on the 
given protein. You can: 

> - Select the phosphosite from the drop-down list of amino acids
> - Type the 3 letter abbreviation of the amino acid and then select the specific one
> - Customise the drop-down list to only contain either SER, THR and/or TYR
> - Type the sequence number of the amino acid you want to select

---

### 2.3 Radius of Interest  

Use the slider to select the radius size from 0Å to 30Å. The graph will automatically update with the new radius value. 

---

### 2.4 Residue Solvent Accessibility Threshold

Use the slider to select Residue Solvent Accessibility Threshold. A larger value (closer to 1)  
means that more of the residue needs to be exposed to the surface to be considered as being 'on the surface'.  
So a smaller value (closer to 0) means that more amino acids are included in the structural motf. 


***

## 3. Viewing Structural Motifs

The main tab is the Visualisation tab which displays the chosen phosphosite on the protein as a simple  
asteroid plot view and also a more abstracted adjacency matrix view. 

---

### 3.1 Asteroid Plot  

The asteroid plot shows all amino acids that are within k hops of the phosphosite. 
Circles that are closer to  
the centre phosphosite are closer in space. Larger circles have a larger Relative Solvent Accessibility, i.e., they  
"stick out" of the protein more. The colour of the amino acid represents its hydrophobicity.

---

### 2.2 Adjacency Matrix

The adjacency matrix is a visualisation of the motif in terms of the distance between the amino acids in the  
protein at a specified radius. 

&nbsp

The **top slider** allows you to choose the radius used for the area of the protein that you wish to include in  
the structural motif. Just slide the circle on the bar from left to right to select any radius between  
0 and 30 angstrom.

&nbsp

The **bottom slider** allows you to adjust the surface accessibility threshold between 0 and 1. This means that  
you can customise which amino acids are considered as being on the surface of the protein.

&nbsp

You can choose to order the amino acids on the axis of the matrix by either **hydrophobicity** or **sequence  
position**. Just use the drop-down menu on the right of the matrix to select the order you wish and the  
matrix will automatically adjust in response.

---

### 3.3 Interactivity

#### Graph Features

If you hover your mouse inside the perimeter of any graph you will be able to see a menu pop up in the  
top right-hand corner. These features are produced by Plotly and they include (moving from left to right):

> (1) Downloading the plot as a png  
> (2) If zoom is selected, then dragging your mouse across the graph means it  zooms in  
> (3) If pan is selected, then dragging your mouse across the graph allows you  to traverse it  
> (4) Zoom in  
> (5) Zoom out  
> (6) Autoscales the graph to its complete size  
> (7) Resets the axes back to the original layout   
> (8) Link to the Plotly website  


&nbsp

#### Grayscale

There is a drop-down option which allows you to customise whether the graphs are in colour or in grayscale  
and you can change it at any point while using the tool. 

***

## 4. Comparing Structural Motifs

The second tab is the Compare Motifs tab which allows you to select the number of structural motifs that  
you want to compare and they will be displayed in a grid format. 


***

## 5. Proteome View

---

### 5.1 Preprocessing

The Proteome View tab allows you to see pre-computed clusterings of structural motifs.  Currently, two sets of proteins  
are supported: Human (known kinases), and Yeast (all). The software is able to handle an arbitrary number of different proteomes,  
but preprocessing must first be completed using the CLI.  

Data from the Phospho.ELM database was used with a list of entries containing  
a. protein accession ID  
b. phosphosite residue ID (sequence position, amino acid type)  
c. kinase label (if known)   

&nbsp

AlphaFold structures were then retrieved and subgraphs constructed, centred around each entry's phosphosite.  
Each subgraph represented a structural motif. 
A Graph Neural Network approach was used to perform unsupervised  
representation learning of the structural motifs in a low-dimension space.  
Embeddings of 80 dimensions were obtained for each structural motif. 


***

> Latest update: 2/8/2022
''')




'''
## Contents:

1. [Input](<h2>1. Input</h2>)
+ 1.1 [Protein](###-1.1-protein)
+ 1.2 [Phosphosite](-1.2-phosphosite)
+ 1.3 [Radius of Interest](-1.2-phosphosite)
2. [Viewing Structural Motifs](##-2.-Viewing-Structural-Motifs)
+ 2.1 [Simple Graph](###-2.1-Simple-Graph)
+ 2.2 [Adjacency Matrix](##-1.-Input)
+ 2.3 [3D Interactive View](##-1.-Input)
3. [Comparing Structural Motifs](##-3.-Comparing-Structural-Motifs)
+ 3.1 [Two Motifs](##-1.-Input)
+ 3.2 [Three or More Motifs](##-1.-Input)
4. [Page Elements](##-1.-Input)
+ 4.1 [Side Panel](##-1.-Input)
+ 4.2 [Different Tabs](##-1.-Input)
5. [Saving](##-1.-Input)

'''

