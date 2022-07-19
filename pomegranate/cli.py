"""Command line interface for pomegranate."""

import pathlib
import click as c

import app
from validate import get_database

@c.command
@c.option(
    "-p",
    "--protein-id",
    help="The accession ID for the protein of interest",
    #default="P0DP23", # Calmodulin 
    #default="4hhb", # PDB accession, NOT UniProt
)
@c.option(
    "--input",
    help="The path to a local .PDB file",
    type=c.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)
def main(protein_id, input):
    
    """Run web application"""
    if not protein_id or input:
        app.run(protein_id)
        return 
    

    protein_id = protein_id.upper()
    if protein_id:
        db = dict(
            uniprot="UniProt",
            pdb="Protein Data Bank",
            unknown="unknown",                
        )
        db_name = get_database(protein_id)
        c.echo(f"Protein {protein_id} from {db[db_name]}")
        
    
if __name__ == "__main__":
    main()