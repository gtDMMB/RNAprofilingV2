
import data
import networkx as nx
import numpy as np
import itertools
from hellinger import hellinger_numpy
from structure_dataframe import StructureDataframe

def rev_shift(helix):
    return (helix[0] - 1, helix[1] - 1, helix[2])

def generate_leaf_radial_diagrams(folder, G, helix_structures, helix_class_labels, sequence):
    
    import draw
    from collections import Counter

    reversed_label_map = data.flip_dict(helix_class_labels)

    for node, node_data in G.nodes(data=True):
        if G.out_degree(node) != 0:
            continue

        print("Radial diagram ", node)
        
        structure_idxs = node_data["structure_idxs"]
        node_structures = [tuple(helix_structures[idx]) for idx in structure_idxs]

        if len(node_structures) == 0:
            print(node_data["count"])
            continue

        structure_counts = Counter(node_structures)
        most_common_struct = max(node_structures, key = structure_counts.get)

        filename = folder + "radial_diagram_{}.svg".format(node)

        draw.plot_radial_diagram(
            most_common_struct, 
            sequence, 
            filename,
            node_data["leaf_label"])

def generate_node_arc_diagrams(folder, G, helix_structures, helix_class_labels, sequence):

    import draw

    reversed_label_map = data.flip_dict(helix_class_labels)

    max_diameter = 0
    for helix_class in helix_class_labels:
        i,j,k = helix_class
        if (j - i) > max_diameter:
            max_diameter = (j - i)

    empty_filename = folder + "arc_diagram_default.svg"

    figure_ratio = max_diameter / len(sequence) * 0.6
    draw.generate_region_arc_diagram(
        len(sequence),
        [],
        [],
        [],
        [],
        figure_ratio,
        empty_filename)

    for node, node_data in G.nodes(data=True):

        print("arc diagram " + str(node))
        
        if G.nodes[node]["type"] != "root":
            parent_node = next(G.predecessors(node))
        else:
            parent_node = None

        parent_edges = list(G.predecessors(node)) + [node]

        path_nodes = nx.ancestors(G, node)
        
        #path_nodes.add(node)

        parent_implications = sum((G.edges[edge]["decision"] for edge in G.subgraph(parent_edges).edges 
            if "decision" in G.edges[edge]), [])
        implications = sum((G.edges[edge]["decision"] for edge in G.subgraph(path_nodes).edges), [])
        t_implications = [feat for feat, present in implications if present]
        f_implications = [feat for feat, present in implications if not present]

        negative_features, important_features, indep_features = [], [], []
        if G.nodes[node]["type"] not in ["leaf", "root"] or (parent_node is not None and G.nodes[parent_node]["type"] != "contingency"):
            negative_features = [feat for feat, present in parent_implications if not present]
            important_features = [feat for feat, present in parent_implications if present]

        elif G.nodes[node]["type"] == "leaf":
            indep_features = sum((list(decision[0]) + list(decision[1])
                for decision in parent_implications), [])

        keep_helices = set()
        for feature in negative_features + important_features + indep_features + t_implications:
            keep_helices.update(
                [rev_shift(helix) for helix in reversed_label_map[feature]])

        important_helices = set()
        for feature in important_features + negative_features:
            important_helices.update(
                [rev_shift(helix) for helix in reversed_label_map[feature]])

        indep_helices = []
        for feature in indep_features:
            indep_helices.append(
                [rev_shift(helix) for helix in reversed_label_map[feature]])

        negative_helices = set()
        for feature in negative_features:
            negative_helices.update(
                [rev_shift(helix) for helix in reversed_label_map[feature]])

        filename = folder + "arc_diagram_{}.svg".format(node)

        label = None
        if "leaf_label" in node_data:
            label = node_data["leaf_label"]

        draw.generate_region_arc_diagram(
            len(sequence),
            keep_helices,
            indep_helices,
            important_helices,
            negative_helices,
            figure_ratio,
            filename,
            label,
            helix_class_labels)

    
def get_decision_text(decision_list):
    decision_text_list = []
    for decision in decision_list:
        if len(decision[1]) == 0:
            decision_text_list.append(", ".join(decision[0]))
            continue

        decision_text_list.append(", ".join(decision[0]) + " || " + ", ".join(decision[1]))

    return ", ".join("(" + text + ")" for text in decision_text_list)


def prepare_agraph_attrs_new(G):
    for edge in G.edges(data = True):
        edge_data = edge[2]
        edge_data["id"] = str(edge[0]) + "_" + str(edge[1])
        edge_data["color"] = "gray18"
        if "decision" in edge_data:
            if "type" in G.nodes[edge[0]] and G.nodes[edge[0]]["type"] == "contingency":
                edge_data["label"] = get_decision_text(edge_data["decision"])
                edge_data["style"] = "dashed"
            else:
                feature_list = []
                present_features = ",".join([feature for feature, present in edge_data["decision"]
                                    if present])
                absent_features =  ",".join(["¬" + feature for feature, present in edge_data["decision"]
                                    if not present])
                
                sep = ","
                if len(present_features) + len(absent_features) > 4:
                    sep = ",\l"
                if len(present_features) > 0 and len(absent_features) > 0:
                    edge_data["label"] = present_features + sep + absent_features
                else:
                    edge_data["label"] = present_features + absent_features

    leaf_idx = 0
    from ThirdParty import roman

    for node in G.nodes(data = True):
        node_data = node[1]

        node_data["label"] = ""
        node_data["id"] = str(node[0])

        label_rows = []
        if "count_old" in node_data:
            label_rows.append(str(node_data["count_old"]))
        elif "count" in node_data:
            label_rows.append(str(node_data["count"]))

        if "type" in node_data and node_data["type"] == "contingency":
            node_data["style"] = "dashed"

        if G.out_degree(node[0]) == 0:
            leaf_idx += 1
            roman_numeral = roman.toRoman(leaf_idx)
            if "count_old" in node_data:
                label_rows.append(roman_numeral + " - " + str(node_data["count"]))
            else:
                label_rows.append(roman_numeral)
            node_data["shape"] = "box"
            node_data["leaf_label"] = roman_numeral

        node_data["label"] = "\n".join(label_rows)

def get_coverage_count(tree, count="count"):
    total = 0
    for node in tree.nodes:
        if tree.out_degree(node) == 0:
            total += tree.nodes[node][count]

    return total

def save_stem_legend_data(
        filename, 
        reversed_feature_dict, 
        feature_counts, 
        helix_counts, 
        helix_class_labels):
    import json

    sorted_keys = sorted(reversed_feature_dict.keys(), key = lambda x: -feature_counts[x])

    grouped_classes = []

    data_table = []
    for key in sorted_keys:
        value = reversed_feature_dict[key]
        region = data.Find_Stem_Region(value)
        label_list = list(sorted(helix_class_labels[helix_class] for helix_class in value))

        data_table.append({
            "Feature": key,
            "Region": str(region),
            "Frequency": feature_counts[key],
            "Helix Classes": ", ".join(str(lab) for lab in label_list)})

        grouped_classes += value

    class_data_table = []
    for helix_class in sorted(grouped_classes, key = lambda x: -helix_counts[x]):
        class_data_table.append({
            "Helix Class": helix_class_labels[helix_class],
            "(i, j, k)": str(helix_class),
            "Frequency": helix_counts[helix_class]})

    json_data = json.dumps(data_table)
    json_class_data = json.dumps(class_data_table)
    with open(filename, "w") as f:

        f.write("var legendJSON = `\n")
        f.write(json_data)
        f.write("\n`;\n")
        f.write("var additionalLegendJSON = `\n")
        f.write(json_class_data)
        f.write("\n`;")

def save_leaf_data(filename, G):
    import json

    file_data = []
    for node in G.nodes():
        if G.out_degree(node) != 0:
            continue

        label = G.nodes[node]["leaf_label"]
        ancestor_list = sorted(list(nx.ancestors(G, node)))

        file_data.append({
            "id":node,
            "label":label,
            "ancestors":ancestor_list})

    json_data =json.dumps(file_data)
    with open(filename, "w") as f:
        f.write("var leafJSON = `\n")
        f.write(json_data)
        f.write("\n`;")

def save_indep_node_data(filename, G, feature_df, label_structures, label_dict, reversed_stem_label_dict):

    import json
    import math
    import itertools

    table_dict = {}
    for node in G.nodes():

        if G.out_degree(node) != 0:
            continue

        parent_node = next(G.predecessors(node))
        if G.nodes[parent_node]["type"] != "contingency":
            continue
        print("Contingency: ", node)

        decision_count = len(G.nodes[parent_node]["decision"])
        current_df = feature_df.get_original_array()[G.nodes[node]["structure_idxs"]]
        column_idx_dict = feature_df.column_idx_dict

        transformed_array = np.zeros((current_df.shape[0], decision_count),dtype=bool)

        for idx, decision in enumerate(G.nodes[parent_node]["decision"]):
            left_implication = get_implications(decision[0], decision[1])

            structure_on_left = np.all(
                [current_df[:,column_idx_dict[feat]] == present for feat, present in left_implication],axis=0)
            transformed_array[:, idx] = structure_on_left

        values, counts = np.unique(transformed_array, return_counts=True, axis=0)
        if len(values.shape) == 1:
            np.expand_dims(values, -1)
        freq_dict = {tuple(value):count/transformed_array.shape[0] for value, count in zip(values, counts)}

        print(values)
        print(counts)
        print(freq_dict)

        row_decisions, col_decisions = [], []
        if decision_count <= 4:
            col_decisions = G.nodes[parent_node]["decision"][:math.floor(decision_count / 2)]
            row_decisions = G.nodes[parent_node]["decision"][math.floor(decision_count / 2):]
        else:
            col_decisions = G.nodes[parent_node]["decision"][:2]
            row_decisions = G.nodes[parent_node]["decision"][2:]

        crosstab_list = [["" for _ in 
                            range(2**len(col_decisions) + len(row_decisions))]
                            for _ in 
                            range(2**len(row_decisions) + len(col_decisions))]

        row_labels = list(itertools.product(*([[True, False]] * len(row_decisions))))
        col_labels = list(itertools.product(*([[True, False]] * len(col_decisions))))

        print(crosstab_list)

        for idx,row_decision_present in enumerate(row_labels):
            for row_decision_idx, present in enumerate(row_decision_present):
                print(row_decisions[row_decision_idx][int(present)])
                label = ",".join(row_decisions[row_decision_idx][int(not present)])
                if label == "":
                    label = ",".join("¬" + feat for feat in row_decisions[row_decision_idx][int(present)])
                print(label)
                crosstab_list[idx + len(col_decisions)][row_decision_idx] = label 

        for idx,col_decision_present in enumerate(col_labels):
            for col_decision_idx, present in enumerate(col_decision_present):
                label = ",".join(col_decisions[col_decision_idx][int(not present)])
                if label == "":
                    label = ",".join("¬" + feat for feat in col_decisions[col_decision_idx][int(present)])
                crosstab_list[col_decision_idx][idx + len(row_decisions)] = label 
            
        for (row_idx,row_decision_present),(col_idx,col_decision_present) in itertools.product(
                enumerate(row_labels), enumerate(col_labels)):
            structure = tuple(col_decision_present + row_decision_present)
            frequency = 0
            if structure in freq_dict:
                frequency = round(freq_dict[structure],3)
            crosstab_list[row_idx + len(col_decisions)][col_idx + len(row_decisions)] = frequency

        table_dict[node] = crosstab_list

        print(crosstab_list)

    json_data = json.dumps(table_dict, indent=2)
    with open(filename, "w") as f:
        f.write("var indepNodeJSON = `\n")
        f.write(json_data)
        f.write("\n`;")
        f.write("\n")

def save_node_sample_indices(filename, G, helix_structures):
    import json
    from collections import Counter

    index_data = {}
    for node, node_data in G.nodes(data=True):

        structure_idxs = node_data["structure_idxs"]
        node_structure_pairs = [(tuple(helix_structures[idx]),int(idx)) for idx in structure_idxs]
        node_structures = [tuple(helix_structures[idx]) for idx in structure_idxs]

        if len(structure_idxs) == 0:
            continue

        structure_counts = Counter([elem for elem in node_structures])
        most_common_struct = max(node_structure_pairs, key = lambda elem: structure_counts[elem[0]])

        row_dict = {
            "mostCommonIndex":most_common_struct[1],
            "allIndices":[int(idx) for idx in structure_idxs]}
        index_data[node] = row_dict

    json_data = json.dumps(index_data, indent=2)
    with open(filename, "w") as f:
        f.write("var nodeSampleIndicesJSON = `\n")
        f.write(json_data)
        f.write("\n`;")

def get_implications(true_features, false_features):
    decision_implications = [(feature, True) for feature in true_features] + \
                            [(feature, False) for feature in false_features]
    return decision_implications

from itertools import product
def get_implication_set(decision_set):

    implication_pairs = []
    for left_features, right_features in decision_set:
        implication_pairs.append(
            (get_implications(left_features, right_features), 
             get_implications(right_features, left_features)))
    
    implication_list = []
    for implications in product(*implication_pairs):
        implication_list.append(tuple(sorted(sum(implications,[]))))

    #print("Formed implication list")
    #print(decision_set, implication_list)
    return set(implication_list)

def merge_contingency_leaves(tree, *args):
    contingency_nodes = [node for node in tree.nodes if tree.nodes[node]["type"] == "contingency"]

    for contingency_node in contingency_nodes:
        
        children = list(tree.successors(contingency_node))
        
        new_arg_values = {}

        for child in children:
            for arg in args:
                if arg not in new_arg_values:
                    new_arg_values[arg] = tree.nodes[child][arg]
                else:
                    new_arg_values[arg] += tree.nodes[child][arg]

        tree.remove_nodes_from(children[1:])

        tree.edges[contingency_node, children[0]]["decision"] = tree.nodes[contingency_node]["decision"]

        for arg in args:
            tree.nodes[children[0]][arg] = new_arg_values[arg]

def _build_tree_new_recr(tree, current_node, current_samples, decision_history_pairs, leaf_cutoff=25):

    forced_cutoff = leaf_cutoff

    if len(current_samples.columns) == 0:

        leaf_node = tree.order()
        tree.add_node(leaf_node)
        tree.add_edge(current_node, leaf_node)

        tree.nodes[leaf_node]["type"] = "leaf"

        return

    current_count = np.sum(current_samples.counts)

    if current_count < forced_cutoff:
        return

    has_forced_feature = False
    forced_decision = None

    #cluster features together to form decisions
    decision_dict = find_decision_dict(current_samples, leaf_cutoff)

    #choose a decision

    score_tuple_list = []
    
    for label, (left, right) in decision_dict.items():

        features = left + right

        left_decision_implications = get_implications(left, right)
        right_decision_implications = get_implications(right, left)
        
        left_samples = current_samples.subset(left_decision_implications)
        right_samples = current_samples.subset(right_decision_implications)

        left_counts = left_samples.counts
        right_counts = right_samples.counts

        score = hellinger_numpy(left_samples, left_counts / np.sum(left_counts), right_samples, right_counts / np.sum(right_counts))

        #check for forced node
        if np.all(left_counts < forced_cutoff) or np.all(right_counts < forced_cutoff):
            if not np.all(left_counts < forced_cutoff):
                score = score + len(current_samples.columns) + sum(left_counts)
            if not np.all(right_counts < forced_cutoff):
                score = score + len(current_samples.columns) + sum(right_counts)

        score_tuple_list.append((score, tuple(sorted(features)), label))

    score_tuple_list.sort()

    best_label = score_tuple_list[-1][-1]
    best_decision = decision_dict[best_label]

    features = sum(best_decision, [])

    left_decision_implications = get_implications(best_decision[0], best_decision[1])
    right_decision_implications = get_implications(best_decision[1], best_decision[0])

    #construct new decision histories for two sides of decision

    left_samples = current_samples.subset(left_decision_implications)
    right_samples = current_samples.subset(right_decision_implications)

    left_counts = left_samples.counts
    right_counts = right_samples.counts

    left_decision_history = decision_history_pairs + left_decision_implications
    right_decision_history = decision_history_pairs + right_decision_implications

    #insert node and edge and recurse for left side of decision
    if left_samples.shape[0] > 0 and not np.all(left_counts < forced_cutoff):
        left_node = tree.order()
        tree.add_node(left_node)
        tree.add_edge(current_node, left_node)
        
        tree.edges[current_node, left_node]["decision"] = left_decision_implications
        tree.nodes[left_node]["count"] = np.sum(left_samples.counts)
        tree.nodes[left_node]["structure_idxs"] = left_samples.index
        tree.nodes[left_node]["type"] = "decision"

        _build_tree_new_recr(tree, left_node, left_samples, left_decision_history, leaf_cutoff)

    #insert node and edge and recurse for right side of decision
    if right_samples.shape[0] > 0 and not np.all(right_counts < forced_cutoff):
        right_node = tree.order()
        tree.add_node(right_node)
        tree.add_edge(current_node, right_node)
        
        tree.edges[current_node, right_node]["decision"] = right_decision_implications
        tree.nodes[right_node]["count"] = np.sum(right_samples.counts)
        tree.nodes[right_node]["structure_idxs"] = right_samples.index
        tree.nodes[right_node]["type"] = "decision"

        _build_tree_new_recr(tree, right_node, right_samples, right_decision_history, leaf_cutoff)

    return

def build_tree_new(samples_df, leaf_cutoff=25):

    tree = nx.DiGraph()

    root = 0;
    tree.add_node(root)

    tree.nodes[root]["type"] = "root"
    tree.nodes[root]["count"] = np.sum(samples_df.counts)
    tree.nodes[root]["structure_idxs"] = samples_df.index

    _build_tree_new_recr(tree, root, samples_df, [], leaf_cutoff)

    return tree

import networkx as nx

def find_decision_dict(sample_df, min_node_freq):

    if len(sample_df.columns) == 1:
        decision_dict = {0:(list(sample_df.columns),[])}
        return decision_dict
    
    paired_sample_dict = {}
    for feature in sample_df.columns:
        left = [feature]
        right = []

        left_implications = get_implications(left, right)
        right_implications = get_implications(right, left)

        left_samples = sample_df.subset(left_implications, frequency_cutoff = min_node_freq)
        right_samples = sample_df.subset(right_implications, frequency_cutoff = min_node_freq)

        left_smaller = left_samples < right_samples
        if left_smaller:
            key = (left_samples, right_samples)
        else:
            key = (right_samples, left_samples)

        if key in paired_sample_dict:
            paired_sample_dict[key].append((left_smaller, feature))
        else:
            paired_sample_dict[key] = [(left_smaller, feature)]
    
    #print("DICT:")
    decision_dict = {}
    for idx, key in enumerate(sorted(paired_sample_dict.keys())):
        #print(key[0].array,key[1].array)
        #print(paired_sample_dict[key])
        if len(paired_sample_dict[key]) == 1:
            decision_dict[idx] = ([paired_sample_dict[key][0][1]],[])
        else:
            decision_dict[idx] = \
                (sorted(feature for on_left, feature in paired_sample_dict[key] if on_left),
                sorted(feature for on_left, feature in paired_sample_dict[key] if not on_left))
    
    return decision_dict
        
def fuzz_structure(label_structure, label_bp_dict, label_count_cutoffs, bp_structure):
    
    available_bps = set(bp_structure) - set().union(*(label_bp_dict[label] for label in label_structure))
    available_keys = label_bp_dict.keys() - set(label_structure)

    new_keys = []
    for key in available_keys:
        key_bps = label_bp_dict[key]
        intersection_bps = available_bps.intersection(key_bps)

        if len(intersection_bps) >= label_count_cutoffs[key]:
            new_keys.append(key)

    res = tuple(list(label_structure) + new_keys)

    return tuple(sorted(res))

def average_region_bps(label_structures, label_bp_dict, bp_structures):

    res = {}
    for key, bps in label_bp_dict.items():
        count_list = [len(set(structure_bps).intersection(bps))
            for labels, structure_bps in zip(label_structures, bp_structures) 
            if key in labels]

        res[key] = sum(count_list) / len(count_list)

    return res

def augment_tree_counts_recr(node, tree, structures, augment_count_label, augment_idx_label):

    tree.nodes[node][augment_count_label] = np.sum(structures.counts)
    if augment_idx_label is not None:
        tree.nodes[node][augment_idx_label] = structures.index

    for child in tree.successors(node):
        
        if "decision" not in tree.edges[node, child]:
            augment_tree_counts_recr(child, tree, structures, augment_count_label, augment_idx_label)
            continue

        child_implications = tree.edges[node, child]["decision"]
        child_samples = structures.subset(child_implications)

        augment_tree_counts_recr(child, tree, child_samples, augment_count_label, augment_idx_label)

def get_edge_decision(tree, edge):
    implications = tree.edges[edge]["decision"]
    decision = [tuple([
        tuple(sorted(feature for feature, present in implications if present)),
        tuple(sorted(feature for feature, present in implications if not present))
    ])]
    return decision

def label_binary_decisions(tree):
    
    for node in tree.nodes:
        if tree.out_degree(node) not in [1, 2]:
            tree.nodes[node]["decision"] = []
            continue
        left_edge = next(iter(tree.out_edges(node)))

        if "decision" not in tree.edges[left_edge]:
            tree.nodes[node]["decision"] = []
            continue

        decision = get_edge_decision(tree, left_edge)

        tree.nodes[node]["decision"] = decision

def augment_tree_counts(tree, structures, augment_count_label, augment_idx_label=None):
    for n, d in tree.in_degree():
        if d == 0:
            root = n
            break

    augment_tree_counts_recr(root, tree, structures, augment_count_label, augment_idx_label)

from itertools import combinations
def overlapping_features(label_bp_dict, all_present_basepairs = None):
    res = []
    for key_a, key_b in combinations(label_bp_dict.keys(), 2):
        bps_a = label_bp_dict[key_a]
        bps_b = label_bp_dict[key_b]

        overlap = bps_a.intersection(bps_b)
        basepair_union = bps_a | bps_b

        if all_present_basepairs is not None:
            overlap = overlap.intersection(all_present_basepairs)
            basepair_union |= all_present_basepairs
        if len(overlap) > 0:
            res.append((key_a, key_b, len(overlap) / len(basepair_union)))

    return res

def get_basepair_coverage(structure_list, basepair_structure_list, feature_to_basepair_dict):
    total_covered_basepairs = 0
    for structure, basepairs in zip(structure_list, basepair_structure_list):
        structure_basepairs = itertools.chain.from_iterable(feature_to_basepair_dict[feature] for feature in structure)
        covered_basepairs = set(structure_basepairs).intersection(basepairs)
        total_covered_basepairs += len(covered_basepairs)

    return total_covered_basepairs

def _merge_complete_binary_trees_new_recr(node, tree, distance_dict, fuzzy_cutoff=0.74):
    
    distance_to_leaves = None
    present_leaf_count = 0
    min_leaf_idx = None
    for child in nx.descendants(tree, node):
        if "type" not in tree.nodes[child]:
            continue

        if tree.nodes[child]["type"] != "leaf":
            continue

        if distance_to_leaves is None:
            distance_to_leaves = distance_dict[node][child]

        if distance_to_leaves != distance_dict[node][child]:
            present_leaf_count = 0
            break

        if min_leaf_idx is None:
            min_leaf_idx = child

        if child < min_leaf_idx:
            min_leaf_idx = child

        present_leaf_count += 1

    if (distance_to_leaves is not None 
            and distance_to_leaves > 1
            and present_leaf_count >= (2 ** (distance_to_leaves - 1)) * fuzzy_cutoff):

        #print("starting node merging")

        decision_set = set()
        for edge in nx.dfs_edges(tree, node):
            if "decision" in tree.edges[edge]:

                decision_t = tuple(sorted(feature for feature, present in tree.edges[edge]["decision"] if present))
                decision_f = tuple(sorted(feature for feature, present in tree.edges[edge]["decision"] if not present))
                decision = tuple(sorted([decision_t, decision_f]))
                decision_set.add(decision)
            
        flattened_decisions = sum((list(decision[0]) + list(decision[1]) for decision in decision_set), [])
        feature_keys, feature_counts = np.unique(flattened_decisions, return_counts=True)

        if np.all(feature_counts == 1):

            leaves = [node for node in nx.descendants(tree, node) if tree.nodes[node]["type"] == "leaf"]

            #print("handling leaves")

            leaf_implications = {}
            for leaf in leaves:
                #print("Leaf: ", leaf)
                path_nodes = nx.ancestors(tree, leaf) & nx.descendants(tree, node)
                path_nodes.add(node)
                path_nodes.add(leaf)

                if path_nodes is None:
                    leaf_implications[leaf] = []
                    continue

                edges = tree.subgraph(path_nodes).edges

                #print(path_nodes)

                implications = tuple(sorted(sum(
                    (tree.edges[edge]["decision"] for edge in edges 
                    if "decision" in tree.edges[edge]),[])))

                leaf_implications[leaf] = implications


            #print(leaf_implications)

            tree.remove_nodes_from(node for node in nx.descendants(tree, node) 
                                    if tree.nodes[node]["type"] != "leaf")

            tree.nodes[node]["type"] = "contingency"
            tree.nodes[node]["decision"] = list(decision if len(decision[0]) > 0 else tuple(reversed(decision)) 
                                                    for decision in decision_set)

            for leaf in leaves:
                tree.add_edge(node, leaf)
                tree.edges[node, leaf]["decision"] = leaf_implications[leaf]

            missing_implications = get_implication_set(decision_set) - set(leaf_implications.values())

            for idx, implication in enumerate(missing_implications):
                new_leaf = "l_" + str(node) + "_" + str(idx + 1)
                tree.add_node(new_leaf)
                tree.add_edge(node, new_leaf)
                tree.edges[node, new_leaf]["decision"] = implication

                tree.nodes[new_leaf]["type"] = "leaf"

            return


    for child in tree.successors(node):
        _merge_complete_binary_trees_new_recr(child, tree, distance_dict)

def merge_complete_binary_trees_new(tree):
    
    for n, d in tree.in_degree():
        if d == 0:
            root = n

    p = nx.shortest_path_length(tree)
    distance_dict = {source:dic for source, dic in p}

    _merge_complete_binary_trees_new_recr(root, tree, distance_dict)

def get_Featured_Helix_Classes(hc_structures):

    hc_counts = data.count_features(hc_structures)
    hc_list = sorted(hc_counts.keys(),
            key = lambda helix_class: -hc_counts[helix_class])

    entropy_cutoff_idx, _ = data.cutoff_objects_by_entropy(
        hc_list, hc_counts)

    featured_hc_list = hc_list[:entropy_cutoff_idx]

    return featured_hc_list

def Fuzz_Stem_Structures(reversed_stem_dict, feat_stem_structures, basepair_structures, return_fuzzy_bp_dict=False):
    stem_region_dict = {key:data.Find_Stem_Region(hc_list) for key, hc_list in reversed_stem_dict.items()}

    region_tol=5
    count_tol=1./3.

    fuzzy_stem_region_dict = {key:(i-region_tol,j+region_tol,k-region_tol,l+region_tol) for key, (i,j,k,l) in stem_region_dict.items()}
    fuzzy_stem_bp_dict = {key:set(data.Region_To_Basepairs(region)) for key, region in fuzzy_stem_region_dict.items()}

    mean_bp_dict = average_region_bps(feat_stem_structures, fuzzy_stem_bp_dict, basepair_structures)
    fuzzy_stem_count_cutoffs = {key:mean_count*count_tol for key, mean_count in mean_bp_dict.items()}

    fuzzy_stem_structures = [fuzz_structure(stem_structure, fuzzy_stem_bp_dict, fuzzy_stem_count_cutoffs, bp_structure)
                                for stem_structure, bp_structure in zip(feat_stem_structures, basepair_structures)]

    if return_fuzzy_bp_dict:
        return fuzzy_stem_structures, fuzzy_stem_bp_dict

    return fuzzy_stem_structures

def build_and_clean_tree(dataframe, min_node_freq, auxilary_dataframe_list = None, auxilary_dataframe_names = None):
    if auxilary_dataframe_list is None:
        auxilary_dataframe_list = []
        auxilary_dataframe_names = []

    tree = build_tree_new(dataframe, min_node_freq)

    label_binary_decisions(tree)

    merge_complete_binary_trees_new(tree)

    augment_tree_counts(tree, dataframe, "count", "structure_idxs")
    for aux_dataframe, aux_name in zip(auxilary_dataframe_list, auxilary_dataframe_names):
        augment_tree_counts(tree, aux_dataframe, aux_name + "_count", aux_name + "_structure_idxs")

    merge_contingency_leaves(tree, "count", "structure_idxs", 
        *[aux_name + "_count" for aux_name in auxilary_dataframe_names],
        *[aux_name + "_structure_idxs" for aux_name in auxilary_dataframe_names])
    remove_empty_forced_edges(tree, ignore_leaves=False)

    return tree

def remove_empty_forced_edges(tree, ignore_leaves=True):

    for n, d in tree.in_degree():
        if d == 0:
            root = n

    contract_edges = []
    for edge in nx.dfs_edges(tree, root):
        if tree.out_degree(edge[0]) != 1:
            continue
        if "decision" in tree.edges[edge] and tree.edges[edge]["decision"] != []:
            continue
        if ignore_leaves and tree.nodes[edge[1]]["type"] == "leaf":
            continue

        contract_edges.append(edge)

    for edge in reversed(sorted(contract_edges)):
        nx.contracted_nodes(tree, edge[0], edge[1], copy=False, self_loops=False)

def main():
    import sys
    import os
    import itertools

    import data
    import Enumerate

    if len(sys.argv) < 2:
        print("Usage: python3 tree.py sequence_file.fasta [index]")
        exit()

    sequence_file = sys.argv[1]

    repeat_idxs = [1]
    if len(sys.argv) > 2:
        repeat_idxs = [sys.argv[2]]

    data.Init_RNA_Seed()

    data_dict = data.load_sample_sequence(
        sequence_file, 
        repeat_idxs=repeat_idxs,
        structure_count=1000,
        cache_folder="cached_structures")

    sequence_name = data_dict["name"]
    sequence = data_dict["sequence"]

    outfile_root = "decision_tree_" + sequence_name

    basepair_structures = next(iter(data_dict["structures"].values()))[0]
    sample_file = next(iter(data_dict["sample_files"].values()))

    helix_structures = [data.Basepairs_To_Helices(structure) 
                        for structure in basepair_structures]
    helix_class_structures = [data.Helices_To_Helix_Classes(structure, sequence) 
                        for structure in helix_structures]

    featured_classes = get_Featured_Helix_Classes(helix_class_structures)

    helix_class_counts = data.count_features(helix_class_structures)
    helix_classes = sorted(helix_class_counts.keys(),
            key = lambda helix_class: -helix_class_counts[helix_class])

    #helix_labels = Enumerate.generate_helix_class_labels(sequence, min_k=1, hairpin_length=3)
    helix_labels = {helix_class:(idx + 1) for idx, helix_class in enumerate(helix_classes)}
    reversed_label_dict = data.flip_dict(helix_labels)

    featured_labels = {hc:helix_labels[hc] 
            for hc in featured_classes}

    stem_dict = data.Find_Stems(featured_classes, featured_labels)
    reversed_stem_dict = data.flip_dict(stem_dict)

    helix_class_structures = [data.Helices_To_Helix_Classes(structure, sequence) 
                                for structure in helix_structures]
    feat_helix_class_structures = [data.Helix_Classes_To_Profiles(structure, featured_classes) 
                                for structure in helix_class_structures]
    feat_stem_structures = [data.Helix_Classes_To_Stems(structure, stem_dict)
                                for structure in feat_helix_class_structures]

    fuzzy_stem_structures = Fuzz_Stem_Structures(reversed_stem_dict, feat_stem_structures, basepair_structures)

    label_structures = [[helix_labels[helix] for helix in struct] 
                        for struct in helix_class_structures]

    stem_counts = data.count_features(feat_stem_structures)
    fuzzy_stem_counts = data.count_features(fuzzy_stem_structures)

    feat_stem_dataframe = StructureDataframe(feat_stem_structures)
    fuzzy_stem_dataframe = StructureDataframe(fuzzy_stem_structures)

    diff_count = sum(1 for fuzz, strict in zip(fuzzy_stem_structures, feat_stem_structures) if fuzz != strict)
    print("Number changed with fuzzy criterion: ", diff_count)

    fuzzy_stem_dataframe_counts_dict = {tuple(key): value for key, value in zip(fuzzy_stem_dataframe, fuzzy_stem_dataframe.counts)}
    _, min_node_freq = data.cutoff_objects_by_entropy(fuzzy_stem_dataframe_counts_dict.keys(), fuzzy_stem_dataframe_counts_dict)

    print(sorted(fuzzy_stem_dataframe_counts_dict.values()))

    print("Min node freq: ", min_node_freq)

    tree = build_and_clean_tree(fuzzy_stem_dataframe, min_node_freq, [feat_stem_dataframe], ["strict"])

    prepare_agraph_attrs_new(tree)
    tree_agraph = nx.nx_agraph.to_agraph(tree)

    output_folder = "output/"
    arc_diagram_folder = output_folder + "arc_diagram/"
    radial_diagram_folder = output_folder + "radial_diagram/"
    output_data_folder = output_folder + "Data/"

    import os
    import shutil

    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    shutil.copytree("site_code",output_folder)

    if not os.path.exists(arc_diagram_folder):
        os.makedirs(arc_diagram_folder)
    if not os.path.exists(radial_diagram_folder):
        os.makedirs(radial_diagram_folder)
    if not os.path.exists(output_data_folder):
        os.makedirs(output_data_folder)

    tree_agraph.write(output_data_folder + "tree.dot")

    coverage = get_coverage_count(tree)
    coverage_prop = coverage / fuzzy_stem_dataframe.shape[0]

    ##########################################
    generate_node_arc_diagrams(
        arc_diagram_folder,
        tree,
        helix_structures,
        stem_dict,
        sequence)
         
    generate_leaf_radial_diagrams(
        radial_diagram_folder,
        tree,
        helix_structures,
        stem_dict,
        sequence)

    #############################################
    
    print("Fuzzy ", fuzzy_stem_counts)
    print(len(fuzzy_stem_structures))

    save_stem_legend_data(
        output_data_folder + "legendJSON.js", 
        reversed_stem_dict, 
        fuzzy_stem_counts, 
        helix_class_counts, 
        helix_labels)
    save_leaf_data(output_data_folder + "leafJSON.js", tree)
    save_indep_node_data(output_data_folder + "indepNodeJSON.js", tree, fuzzy_stem_dataframe, label_structures, helix_labels, reversed_stem_dict)
    save_node_sample_indices(output_data_folder + "nodeSampleIndicesJSON.js", tree, helix_structures)

    import subprocess
    subprocess.run(["dot","-T","svg","-o", output_data_folder + "tree.svg",output_data_folder + "tree.dot"])

    with open(output_data_folder + "treeSVG.js","w") as f:
        f.write("var treeTXT = `\n")
        with open(output_data_folder + "tree.svg","r") as svg_file:
            svg_text = svg_file.read()
        f.write(svg_text)
        f.write("\n`;")

    with open(output_data_folder + "sequenceName.js","w") as f:
        f.write("var sequenceName = \"" + sequence_name + "\";\n")

    with open(output_data_folder + "sampleGTBOLTZ.js","w") as f:
        f.write("var sampleGTBOLTZ = `\n")
        with open(sample_file,"r") as sample_file:
            next(sample_file)
            sample_text = sample_file.read()
        f.write(sample_text)
        f.write("\n`;")
        
if __name__ == "__main__":
    main()
