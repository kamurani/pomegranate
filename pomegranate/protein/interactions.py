import networkx as nx

from graphein.protein.utils import filter_dataframe
from graphein.protein.edges.distance import compute_distmat, get_interacting_atoms


'''
From graphein.protein.edges.distance
Modified by Naomi Warren
'''
def add_distance_threshold(
    G: nx.Graph, long_interaction_threshold: int, threshold: float = 5.0
):
    """
    Adds edges to any nodes within a given distance of each other.
    Long interaction threshold is used to specify minimum separation in sequence
    to add an edge between networkx nodes within the distance threshold
    :param G: Protein Structure graph to add distance edges to
    :type G: nx.Graph
    :param long_interaction_threshold: minimum distance in sequence for two
        nodes to be connected
    :type long_interaction_threshold: int
    :param threshold: Distance in angstroms, below which two nodes are connected
    :type threshold: float
    :return: Graph with distance-based edges added
    """
    pdb_df = filter_dataframe(
        G.graph["pdb_df"], "node_id", list(G.nodes()), True
    )
    dist_mat = compute_distmat(pdb_df)
    interacting_nodes = get_interacting_atoms(threshold, distmat=dist_mat)
    interacting_nodes = list(zip(interacting_nodes[0], interacting_nodes[1]))

    #log.info(f"Found: {len(interacting_nodes)} distance edges")
    count = 0
    for a1, a2 in interacting_nodes:

        # Don't bother adding self-loops
        if a1 == a2:
            continue

        n1 = pdb_df.at[a1, "node_id"]
        n2 = pdb_df.at[a2, "node_id"]
        n1_chain = pdb_df.at[a1, "chain_id"]
        n2_chain = pdb_df.at[a2, "chain_id"]
        n1_position =pdb_df.at[a1, "residue_number"]
        n2_position = pdb_df.at[a2, "residue_number"]

        condition_1 = n1_chain == n2_chain
        condition_2 = (
            abs(n1_position - n2_position) < long_interaction_threshold
        )

        if not (condition_1 and condition_2):
            count += 1
            if G.has_edge(n1, n2):
                G.edges[n1, n2]["kind"].add("distance_threshold")
            else:
                G.add_edge(n1, n2, kind={"distance_threshold"})
    # log.info(
    #     f"Added {count} distance edges. ({len(list(interacting_nodes)) - count} removed by LIN)"
    # )

    return G