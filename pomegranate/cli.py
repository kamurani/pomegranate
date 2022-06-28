"""Command line interface for pomegranate."""

import argparse

parser = argparse.ArgumentParser(description="""POMEGRANATE
PhOsphosite Motif Explorer - GRAph Network Abstraction Through Embeddings""")

parser.add_argument('database', metavar='database', type=str, nargs=1,
                    help='the database to search for the protein from (AlphaFold or PDB)')
parser.add_argument('PDB_filename', metavar='PDB_filename', type=str,       
                    nargs=1, help='the name of a local PDB file you wish to use')
parser.add_argument('radius', metavar='radius', type=float, nargs=1,
                    help='the radius to be examined around the phosphorylation site')

args = parser.parse_args()