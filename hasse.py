
import networkx as nx
from itertools import combinations

def build_raw_hasse_diagram(profile_list):
    previous_intersection_size = 0
    intersection_profiles = set(profile_list)

    while previous_intersection_size < len(intersection_profiles):
        new_profiles = []
        for prof_a, prof_b in combinations(intersection_profiles, 2):
            intersection = set(prof_a).intersection(set(prof_b))
            new_profiles.append(tuple(sorted(intersection)))

        previous_intersection_size = len(intersection_profiles)
        intersection_profiles |= set(new_profiles)

    G = nx.DiGraph()

    edge_list = []
    for prof_a, prof_b in combinations(intersection_profiles, 2):
        if set(prof_b).issubset(set(prof_a)):
            edge_list.append((prof_b, prof_a))
        if set(prof_a).issubset(set(prof_b)):
            edge_list.append((prof_a, prof_b))

    intersection_profiles = [tuple(sorted(prof)) for prof in intersection_profiles]
    
    G.add_nodes_from(intersection_profiles)
    G.add_edges_from(edge_list)

    G = nx.algorithms.dag.transitive_reduction(G)
    
    return G

