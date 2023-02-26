
def find_bracket(helix_list, helix_labels):
    helix_tree = build_tree(helix_list)
    bracket_string = tree_to_bracket(helix_tree, helix_labels)
    return bracket_string

def find_shape(helix_list):
    helix_tree = build_tree(helix_list)
    helix_tree = simplify_tree(helix_tree)
    shape_string = tree_to_shape(helix_tree)

    return shape_string
    
def get_root_key():
    return (-2, float('inf'), 1)

def insert_helix(helix_tuple, input_key, helix_tree):
    pseudoknot = False

    for test_key in helix_tree[input_key]["children"]:
        if helix_tuple[0] >= test_key[0] and helix_tuple[1] <= test_key[1]:
            return insert_helix(helix_tuple, test_key, helix_tree)
        elif helix_tuple[0] >= test_key[0] and helix_tuple[0] <= test_key[1]:
            pseudoknot = True
        elif helix_tuple[1] >= test_key[0] and helix_tuple[1] <= test_key[1]:
            pseudoknot = True

    helix_tree[input_key]["children"].append(helix_tuple)
    helix_tree[helix_tuple] = {"parent":input_key, "children":[]}

    return pseudoknot

def build_tree(helix_list, check_for_pseudoknot = False):
    
    helix_list = [(x[0], -x[1], x[2]) for x in helix_list]
    helix_list.sort()
    helix_list = [(x[0], -x[1], x[2]) for x in helix_list]

    root_key = get_root_key()
    helix_tree = {}
    helix_tree[root_key] = {"parent":"None", "children":[]}

    pseudoknot = False
    for helix_tuple in helix_list:
        pseudoknot_present = insert_helix(helix_tuple, root_key, helix_tree)
        pseudoknot = (pseudoknot or pseudoknot_present)

    if check_for_pseudoknot:
        return helix_tree, pseudoknot

    return helix_tree

def simplify_tree(helix_tree):

    delete_list = []
    for helix_ijk in helix_tree.keys():
        if len(helix_tree[helix_ijk]["children"]) == 1:
            delete_list.append(helix_ijk)
            
            parent = helix_tree[helix_ijk]["parent"]
            child = helix_tree[helix_ijk]["children"][0]

            if parent != "None":
                helix_tree[parent]["children"].remove(helix_ijk)
                helix_tree[parent]["children"].append(child)

            helix_tree[child]["parent"] = parent

    for helix_ijk in delete_list:
        if helix_ijk in helix_tree:
            del helix_tree[helix_ijk]

    return helix_tree

def subtree_to_shape(helix_tree, current_node):
    output_string = "[" 

    for node in helix_tree[current_node]["children"]:
        output_string += subtree_to_shape(helix_tree, node)

    output_string += "]"
    return output_string

def tree_to_shape(helix_tree):
    start_node = None
    for node in helix_tree:
        if helix_tree[node]["parent"] == "None":
            start_node = node
            break

    output_string = ""
    if start_node != get_root_key():
        output_string += "["


    for node in helix_tree[start_node]["children"]:
        output_string += subtree_to_shape(helix_tree, node)

    if start_node != get_root_key():
        output_string += "]"

    return output_string

def subtree_to_bracket(helix_tree, current_node, parent_label, helix_labels):
    output_string = ""

    label = helix_labels[current_node]
    if label != parent_label and label != "":
        output_string += "[{}".format(label)

    for node in helix_tree[current_node]["children"]:
        output_string += subtree_to_bracket(helix_tree, node, label, helix_labels)

    if label != parent_label and label != "":
        output_string += "]"
    return output_string

def tree_to_bracket(helix_tree, helix_labels):
    start_node = None
    for node in helix_tree:
        if helix_tree[node]["parent"] == "None":
            start_node = node
            break

    label = ""
    if start_node != get_root_key():
        label = helix_labels[start_node]

    output_string = ""
    if start_node != get_root_key():
        output_string += "[{}".format(label) 

    for node in helix_tree[start_node]["children"]:
        output_string += subtree_to_bracket(helix_tree, node, label, helix_labels)

    if start_node != get_root_key():
        output_string += "]"

    return output_string

if __name__ == "__main__":
    
    helix_list = [#(1, 100, 2), (5, 95, 2), 
                (20, 30, 2), (40, 80, 2), 
                (31, 39, 2), (44, 60, 2),
                (65, 75, 2)]

    helix_list = [(1, 100, 2)]

    shape = find_shape(helix_list)
    print(helix_list)
    print(shape) 
    
