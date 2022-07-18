"""Load graphs from list of PDB files"""


from graphein.protein.graphs import construct_graph
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.distance import *	

import os
import urllib as ul

from definitions import STRUCTURE_PATH
from validate import get_database

def load_graphs(
    pdb_path,       # directory containing local pdb files
    pdb_id_list,    # file containing list of psites 
):

    # for each entry: 
    accs = []
    for acc in accs:

        if os.path.exists(pdb_path):
            pass
            #acc 
            #g = 
        else:
            
            acc = "P06738"

            filename = acc + ".pdb"
            url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v2.pdb"

            ul.urlretrieve(url, STRUCTURE_PATH+f"/{filename}")

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