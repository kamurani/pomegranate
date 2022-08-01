import os

# Pomegranate  root 
source_dir = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(source_dir)
STRUCTURE_PATH = os.path.join(ROOT_DIR, 'pomegranate/structures/')

SAVED_GRAPHS_DIR = os.path.join(ROOT_DIR, "graphs")

SAVED_PDB_DIR = os.path.join(ROOT_DIR, 'examples/pdbs/')

STRUCTURE_HUMAN_PATH = os.path.join(STRUCTURE_PATH, 'human')

EMBEDDINGS_FILENAME = "embeddings_output.csv" #'embeddings.csv'
EMBEDDINGS_PATH = os.path.join(os.path.join(ROOT_DIR, 'embeddings'), EMBEDDINGS_FILENAME)

