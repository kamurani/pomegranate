"""Util functions for amino acids."""

from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1


"""
Convert a 3 letter AA code to 1 letter; 
Or itself (if already 1 letter code)
"""
def aa1letter(
    aa: str, 
):
    assert(aa is not None)
    aa = aa.upper()
    if len(aa) == 1:
        return aa
    elif len(aa) == 3:
        try: return aa3to1(aa)
        except: raise ValueError(f"Invalid amino acid '{aa}' given.")
    else: raise ValueError(f"Invalid amino acid '{aa}' given.")

