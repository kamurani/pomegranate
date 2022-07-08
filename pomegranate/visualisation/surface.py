
from pymol import cmd

# Helper functions
def unique_list (res_list):
    npify = np.array(res_list)
    npuniq = np.unique(npify)
    return npuniq.tolist()

# Given a residue and a radius, find all AAs in that region
# Input:  site   (integer -> position in protein sequence)
#         radius (integer -> radius in Angstroms)
# Output: list of residues within [radius] Angstroms of [site]
def res_within (site, radius):
    # Select the desired region
    cmd.select('region', f'(all within {radius} of i. {site})')

    # Iterate over atoms wtihin selection
    all_ress = {'within_r': []}
    cmd.iterate('region', 'within_r.append(resv)', space=all_ress)

    # Remove double-ups
    all_ress['within_r'] = unique_list(all_ress['within_r'])

    return all_ress['within_r']


# Find out the angles between them

# Display 