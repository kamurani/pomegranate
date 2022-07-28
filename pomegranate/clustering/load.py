"""Load graphs from list of PDB files"""


from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

from pathlib import Path





import pandas as pd

import os
import urllib as ul

from definitions import STRUCTURE_PATH
from validate import get_database

# TODO: use Path types instead of strs
def load_graphs(
    pdb_path: str = None,       # directory containing local pdb files
    psite_list: str = None,          # path to file containing list of psites
):

    psite_path = Path(psite_list)
    if not psite_path.is_file():
        raise ValueError(f"No such file {psite_list}")
    
    df = pd.read_csv(psite_path, sep='\t', header=0)

    



    # Default directory
    if pdb_path == None:
        pdb_path = STRUCTURE_PATH + "/yeast"

    if not os.path.isdir(pdb_path):
        raise ValueError(f"Path {pdb_path} is not a directory.")
    
    # for each entry: 
    accs = ["P06738"]
    for acc in accs:

        filename = pdb_path + acc + ".pdb"
        if not os.path.exists(filename):
            
            print(f"Downloading {acc} from AF2...")

            
    
        filename = acc + ".pdb"
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb"

        ul.urlretrieve(url, pdb_path+f"/{filename}")

        # Phosphosite 
        #subgraph = 


    
def main():
    
    

    
    print(STRUCTURE_PATH)
    acc = "P06738"

    
    
    print(get_database(acc))
    return 

    

    load_graphs(
        pdb_path = STRUCTURE_PATH,
    )


if __name__ == "__main__":

    main()