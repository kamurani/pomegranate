import graphein.protein as gp
from graphein.testing import graphs_isomorphic
import networkx.readwrite as nx
import json
# import sys # For cmd line debugging


def g_to_json(g):

    del g.graph["config"] # Remove the config from the graph as it's not easily serialisable

    # Ensure graph data is in JSON format
    for k, v in g.graph.items():
        try:
            g.graph[k] = v.to_json()
        except AttributeError:
            try:
                g.graph[k] = v.tolist()
            except AttributeError:
                continue

    # Ensure node data is in JSON format
    for n, d in g.nodes(data=True):
        for k, v in d.items():
            try:
                d[k] = v.to_json()
            except AttributeError:
                try:
                    d[k] = v.tolist()
                except AttributeError:
                    continue

    # Ensure edge data is in JSON format
    for _, _, d in g.edges(data=True):
        for k, v in d.items():
            try:
                d[k] = v.to_json()
            except AttributeError:
                try:
                    d[k] = v.tolist()
                except AttributeError:
                    try:
                        d[k] = list(v)
                    except AttributeError:
                        continue

    j_graph = nx.json_graph.node_link_data(g)
    return j_graph
    # print(json.dumps(j_graph))
    # # Write the graph to a JSON file
    # with open("test.json", 'w') as f:
    #     tmp = nx.json_graph.node_link_data(g)
    #     f.write(json.dumps(tmp))

    # # Read the graph back in
    # with open("test.json", "r") as f:
    #     json_data = json.load(f)

    # # NB, graph node and edge features are still in JSON form and will need to be converted back into pandas/numpy objects
    # h = nx.json_graph.node_link_graph(json_data)

    # assert graphs_isomorphic(g, h), "Graphs are not isomorphic"

'''
DEBUGGING
'''

# code = sys.argv[1]
# print(code)
#g = gp.construct_graph(pdb_code=code)
# g_to_json(g)
