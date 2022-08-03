from typing import Dict

HYDROPHOBICITY_SCALES: Dict[str, Dict[str, float]] = {
    "kd": { # kdHydrophobicity (a) 
        "ILE": 4.5,
        "VAL": 4.2,
        "LEU": 3.8,
        "PHE": 2.8,
        "CYS": 2.5,
        "MET": 1.9,
        "ALA": 1.8,
        "GLY": -0.4,
        "THR": -0.7,
        "SER": -0.8,
        "TRP": -0.9,
        "TYR": -1.3,
        "PRO": -1.6,
        "HIS": -3.2,
        "GLU": -3.5,
        "GLN": -3.5,
        "ASP": -3.5,
        "ASN": -3.5,
        "LYS": -3.9,
        "ARG": -4.5,
    },
    "ww": { # wwHydrophobicity (b)
        "ILE": 0.31,
        "VAL": -0.07,
        "LEU": 0.56,
        "PHE": 1.13,
        "CYS": 0.24,
        "MET": 0.23,
        "ALA": -0.17,
        "GLY": -0.01,
        "THR": -0.14,
        "SER": -0.13,
        "TRP": 1.85,
        "TYR": 0.94,
        "PRO": -0.45,
        "HIS": -0.96,
        "GLU": -2.02,
        "GLN": -0.58,
        "ASP": -1.23,
        "ASN": -0.42,
        "LYS": -0.99,
        "ARG": -0.81,
    },
    "hh": { # hhHydrophobicity (c)
        "ILE": -0.60,
        "VAL": -0.31,
        "LEU": -0.55,
        "PHE": -0.32,
        "CYS": -0.13,
        "MET": -0.10,
        "ALA": 0.11,
        "GLY": 0.74,
        "THR": 0.52,
        "SER": 0.84,
        "TRP": 0.30,
        "TYR": 0.68,
        "PRO": 2.23,
        "HIS": 2.06,
        "GLU": 2.68,
        "GLN": 2.36,
        "ASP": 3.49,
        "ASN": 2.05,
        "LYS": 2.71,
        "ARG": 2.58,
    },
    "mf": { # mfHydrophobicity (d)
        "ILE": -1.56,
        "VAL": -0.78,
        "LEU": -1.81,
        "PHE": -2.20,
        "CYS": 0.49,
        "MET": -0.76,
        "ALA": 0.0,
        "GLY": 1.72,
        "THR": 1.78,
        "SER": 1.83,
        "TRP": -0.38,
        "TYR": -1.09,
        "PRO": -1.52,
        "HIS": 4.76,
        "GLU": 1.64,
        "GLN": 3.01,
        "ASP": 2.95,
        "ASN": 3.47,
        "LYS": 5.39,
        "ARG": 3.71,
    },
    "tt": { # ttHydrophobicity (e)
       "ILE": 1.97,
        "VAL": 1.46,
        "LEU": 1.82,
        "PHE": 1.98,
        "CYS": -0.30,
        "MET": 1.40,
        "ALA": 0.38,
        "GLY": -0.19,
        "THR": -0.32,
        "SER": -0.53,
        "TRP": 1.53,
        "TYR": 0.49,
        "PRO": -1.44,
        "HIS": -1.44,
        "GLU": -2.90,
        "GLN": -1.84,
        "ASP": -3.27,
        "ASN": -1.62,
        "LYS": -3.46,
        "ARG": -2.57, 
    }
}
"""
Set of (5) dictionaries that map amino acid 3-letter codes to their hydrophobicity. 

The scales included are from Chimera (UCSF) https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html 
and are as follows:

    * kdHydrophobicity
        (a) A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32. 
        https://www.ncbi.nlm.nih.gov/pubmed/7108955 

    * wwHydrophobicity
        (b) Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Wimley WC, White SH. Nat Struct Biol. 1996 Oct;3(10):842-8.
        https://www.ncbi.nlm.nih.gov/pubmed/8836100

    * hhHydrophobicity
        (c) Recognition of transmembrane helices by the endoplasmic reticulum translocon. Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81, supplementary data.
        https://www.ncbi.nlm.nih.gov/pubmed/15674282

        In this scale, negative values indicate greater hydrophobicity. 

    * mfHydrophobicity
        (d)  Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Moon CP, Fleming KG. Proc Natl Acad Sci USA. 2011 Jun 21;108(25):10174-7, supplementary data.
        https://www.ncbi.nlm.nih.gov/pubmed/21606332 

        In this scale, negative values indicate greater hydrophobicity.
    
    * ttHydrophobicity
        (e) An amino acid “transmembrane tendency” scale that approaches the theoretical limit to accuracy for prediction of transmembrane helices: relationship to biological hydrophobicity. Zhao G, London E. Protein Sci. 2006 Aug;15(8):1987-2001.
        https://www.ncbi.nlm.nih.gov/pubmed/16877712 
"""