"""Validate input data"""
import re

def get_database(id=None):

    UNIPROT_ID = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    PDB_ID = "[0-9][A-Z0-9]{3}"

    id = id.upper()

    if re.match(UNIPROT_ID, id):
        return "uniprot"
    if re.match(PDB_ID, id):
        return "pdb"
    else:
        return "unknown"