from dash import dcc

def help_text():
    return dcc.Markdown(''' 

## Need Help?
---
### Welcome to POMEGRANATE!  

POMEGRANATE is an interactive protein structural motif explorer which
allows you to do many things including: 

> - Generating structural motifs of input proteins around a given phosphorylation site
> - Viewing structural motifs in colour and greyscale with interactive information  
> - Comparing two structural motifs 

*(Add more info at each sprint)* 

***

## 1. Input

There are 3 different elements of input required to generate the structural motif: the protein, the specific  
phosphorylation site (phosphosite) on that protein, and the radius of interest.  

*(Add screenshots once the layout is complete)* 

---

### 1.1 Protein 

In order to generate a structural motif of a phosphorylation site, you will need to choose or upload the file  
of the protein that you are interested in.  

&nbsp  

###### Uploading a Protein File
The file either needs to be a PBD or MMCIF file. Just select the "Upload file" button and select the protein file you  
want from your local computer files.  

&nbsp

###### Select from Database with Protein ID
Enter the ID of the protein that you wish to select. The tool will then retrieve the file of the protein you specified  
from the relevant database. 

---

### 1.2 Phosphosite 

There are many options of how to select the phosphosite of interest on the 
given protein. You can: 

> - Select the phosphosite from the drop-down list of amino acids
> - Type the 3 letter abbreviation of the amino acid and then select the specific one
> - Customise the drop-down list to only contain either SER, THR and/or TYR
> - Type the sequence number of the amino acid you want to select

---

### 1.3 Radius of Interest  

Sprint 3

***

## 2. Viewing Structural Motifs

---

### 2.1 Simple Graph  

Asteroid Plot?

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

### 2.3 3D Interactive View

*(describe this more once it is added to the web interface)*  

This 3D visualisation allows you to **save** custom snips of the structure you are observing. 

***

## 3. Comparing Structural Motifs

---

### 3.1 Two Motifs

Select the side-by-side comparison.  
*(need to work out these options)*

---

### 3.2 Three or More Motifs

Grid-view arranged by similarity.

***

## 4. Page Elements

---

### 4.1 Side Panel

*(add once completed)*

----

### 4.2 Different Tabs

The **'Visualisation' tab** is the main screen for the tool where you are able to enter the required input and  
view the structural motif generated from the protein.

&nbsp

The **'Compare Motifs' tab** allows you to view 2 motifs in a side-by-side layout so that you can visually  
compare the similarities and differences. 

&nbsp

The **'Proteome View' tab** ...

***

## 5. Saving

---

Allows you to save screenshots?  

Saving your history as well?

***

> Latest update: Middle of Sprint 2 - 19/7/2022
''')