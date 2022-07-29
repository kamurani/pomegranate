import os

# Pomegranate  root 
source_dir = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(source_dir)
STRUCTURE_PATH = os.path.join(ROOT_DIR, 'structures')
<<<<<<< HEAD
=======

SAVED_GRAPHS_DIR = os.path.join(ROOT_DIR, "graphs")

STRUCTURE_HUMAN_PATH = os.path.join(STRUCTURE_PATH, 'human')

EMBEDDINGS_FILENAME = "embeddings_output.csv" #'embeddings.csv'
EMBEDDINGS_PATH = os.path.join(os.path.join(ROOT_DIR, 'embeddings'), EMBEDDINGS_FILENAME)

>>>>>>> 699e6d4481baa3c2382dc4ff8a01d0654ff7e642
