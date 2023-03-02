
import os
import networkx as nx
from collections import deque

#helices are 0 indexed, helix classes are 1 indexed in this file

from itertools import chain, combinations, product

#modified from itertools docs
def partial_powerset(iterable, max_elements=None):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    if max_elements is None:
        max_elements = len(s)
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(max_elements+1))

def cutoff_objects_by_entropy(object_list, frequency_dict):
    from math import log

    object_list = sorted(object_list, 
            key = lambda thing: -frequency_dict[thing])

    norm = max(frequency_dict.values())
    entropy = 0
    prev_average_entropy = 0

    #print(frequency_dict)
    #print(object_list)

    if norm <= 0:
        return len(object_list), 0

    for idx, thing in enumerate(object_list):
        if frequency_dict[thing] == 1:
            if idx >= 1:
                return idx, frequency_dict[object_list[idx-1]]
            return idx, 1

        fraction = frequency_dict[thing] / norm
        entropy -= fraction * log(fraction)

        if fraction < 1: #equivalent to fraction != 1
            entropy -= (1 - fraction) * log(1 - fraction)

        average_entropy = entropy / (idx + 1)

        #print(average_entropy)

        if (prev_average_entropy < average_entropy or 
                abs(prev_average_entropy - average_entropy) < 0.0000000000000001):
            prev_average_entropy = average_entropy
        else:
            return idx, frequency_dict[object_list[idx-1]]

    return len(object_list), 0

def To_Dot_Bracket(helix_class_list, sequence_length):
    characters = ["."] * sequence_length

    skipped_pairs = []
    for helix_class in sorted(helix_class_list):
        i, j, k = helix_class
        for idx in range(k):
            if characters[i + idx] != "." or\
                characters[j - idx] != ".":

                skipped_pairs.append((i + idx, j - idx))
                continue

            characters[i + idx] = "("
            characters[j - idx] = ")"

    return "".join(characters), set(skipped_pairs)

def Dot_to_BP(dot_string):
    opened_set = set(["(","<","{","["])
    closed_set = set([")",">","}","]"])

    bp_stack = []
    bp_list = []
    for idx, character in enumerate(dot_string):

        if character in opened_set:
            bp_stack.append(idx)

        if character in closed_set:
            matching_idx = bp_stack.pop()
            bp_list.append((matching_idx, idx))

    return bp_list

#seed_set = False
def Init_RNA_Seed(seed=''):
    import RNA

    RNA.init_rand(seed)

    #seed_set = True

def Sample_Sequence(sequence, structure_count=100):
    # load RNAlib library python bindings (associated with the ViennaRNA packages)
    import RNA

    # create model details
    md = RNA.md()
     
    # activate unique multibranch loop decomposition
    md.uniq_ML = 1
     
    # create fold compound object
    fc = RNA.fold_compound(sequence, md)
     
    # compute MFE
    (ss, mfe) = fc.mfe()
     
    # rescale Boltzmann factors according to MFE
    fc.exp_params_rescale(mfe)
     
    # compute partition function to fill DP matrices
    fc.pf()

    structures = fc.pbacktrack(structure_count)

    return structures

def Sample_Sequence_RNAstructure(RNAstructure_location, seed, sequence_file, output_file, structure_count=1000):
        import subprocess
        from pathlib import Path
        import rna_shape

        data_tables = Path(RNAstructure_location) / "data_tables"
        stochastic_program = Path(RNAstructure_location) / "exe/stochastic"
        ct_file = Path("./.tmp.ct")

        subprocess.run([
            stochastic_program.resolve(), 
            "--sequence",
            "--seed", str(seed), 
            "--ensemble", str(structure_count),
            sequence_file, ct_file.resolve()],
            env=dict(os.environ,DATAPATH=data_tables.resolve()))

        structures, sequence_tuple = Read_CT(ct_file.resolve())

        helix_structures = [Basepairs_To_Helices(struct) for struct in structures]
        for idx, helix_struct in enumerate(helix_structures):
            if rna_shape.build_tree(helix_struct, check_for_pseudoknot=True)[1]:
                print("CT file contains a pseudoknot! Pseudoknots are not currently supported. Exiting")
                print("Structure {} was ".format(idx + 1), helix_struct)
                quit()

        if output_file is not None:
            dot_structures = [To_Dot_Bracket(Basepairs_To_Helices(struct), len(sequence_tuple))[0]
                for struct in structures] 

            Write_Dot_Structures(output_file, "".join(sequence_tuple), dot_structures)

        return structures


'''
def Load_Sample_Sequences_In_Folder(folder, repeat_count=1, structure_count=1000, cache_folder="cached_structures"):
    from pathlib import Path
    from itertools import chain

    result_dict = {}

    Path(cache_folder).mkdir(parents=True, exist_ok=True)

    for filename in chain(
            Path(folder).rglob("*.fasta"),
            Path(folder).rglob("*.FASTA")):
        sequence, sequence_name = Read_FASTA(filename)

        result_dict[sequence] = {
            "name": sequence_name,
            "structures": {}}

        for repeat in range(repeat_count):

            cache_file = Path.joinpath(Path(cache_folder), 
                "{}_{}_{}.ct".format(sequence_name, structure_count, repeat))

            if cache_file.is_file():
                structures = Read_CT(cache_file.resolve())
            else:
                structures = Sample_Sequence(sequence, structure_count)
                Write_CT(cache_file.resolve(), sequence, sequence_name, structures)
                structures = [structures, list(sequence)]

            result_dict[sequence]["structures"][repeat] = structures

    return result_dict
'''

def get_sequences_in_folder(folder):
    from pathlib import Path
    from itertools import chain
    
    return list(chain(Path(folder).rglob("*.fasta"),
                    Path(folder).rglob("*.FASTA")))

def get_hash(sequence):
    from hashlib import blake2b

    h = blake2b(digest_size=6)
    h.update(sequence.encode('ascii','ignore'))
    seq_hash = h.hexdigest()

    return seq_hash

def load_sample_sequence(
        sequence_file, 
        repeat_count=1, 
        seed=None,
        structure_count=1000, 
        cache_folder=None,
        RNAstructure_location = None):
    from pathlib import Path

    if cache_folder is not None:
        Path(cache_folder).mkdir(parents=True, exist_ok=True)

    sequence, sequence_name = Read_FASTA(sequence_file)

    if seed is None:
        seed = "NA"

    seq_hash = get_hash(sequence)

    result_dict = {
        "sequence": sequence,
        "name": sequence_name,
        "hash": seq_hash,
        "seed": seed,
        "structures": None,
        "sample_files": None}

    if RNAstructure_location is None:
        sampler = "RNAlib"
    else:
        sampler = "RNAstructure"

    cache_file = Path("")
    if cache_folder is not None:
        cache_file = Path.joinpath(Path(cache_folder), 
            "{}_{}_{}_{}_{}.dot_struct".format(sequence_name, seq_hash, structure_count, seed, sampler))

    if cache_file.is_file():
        dot_structures = Read_Dot_Structures(cache_file.resolve())
        structures = [Dot_to_BP(dot_struct) for dot_struct in dot_structures]
    elif RNAstructure_location is None:
        dot_structures = Sample_Sequence(sequence, structure_count)
        structures = [Dot_to_BP(dot_struct) for dot_struct in dot_structures]

        if cache_folder is not None:
            Write_Dot_Structures(cache_file.resolve(), sequence, dot_structures)
    else:
        structures = Sample_Sequence_RNAstructure(RNAstructure_location, seed, sequence_file, cache_file.resolve(), structure_count)


    result_dict["structures"] = structures

    if cache_folder is not None:
        result_dict["sample_files"] = cache_file

    return result_dict

def Read_FASTA(file_name):

    with open(file_name,'r') as f:
        first = next(f).strip()
        result = ""
        for line in f:
            result += line.strip()

        sequence = result.upper().replace("T","U")
        name = first.replace(">","").split(" ")[-1].replace("/","-")
    
    return sequence, name

def Write_Dot_Structures(file_name, sequence, dot_structures):
    with open(file_name, "w") as f:
        f.write("> {}\n".format(sequence))
        for struct in dot_structures:
            f.write("{}\n".format(struct))

def Read_Dot_Structures(file_name):
    result = []
    with open(file_name, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            result.append(line.strip().split()[0])

    return result

def Write_CT(file_name, sequence, sequence_name, structures):
    with open(file_name, 'w') as f:
        for structure in structures:
            f.write("{:>5}  {}\n".format(len(sequence), sequence_name))

            bp_dict = {}
            for bp in structure:
                bp_dict[bp[0]] = bp[1]
                bp_dict[bp[1]] = bp[0]
            for idx in range(len(sequence)):
                pair_elem = 0
                if idx in bp_dict:
                    pair_elem = bp_dict[idx] + 1

                f.write(" {:>4} {} {:>7} {:>4} {:>4} {:>4}\n".format(
                    idx + 1, 
                    sequence[idx],
                    idx,
                    idx + 2,
                    pair_elem,
                    idx + 1))

# Reads the base pairs from a ct file
# Returns a list of lists of tuples where each tuple is a basepair and
# also returns the sequence of bases
def Read_CT(file_name):
    single_structure = []
    structure_list = []

    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            line_data = [x.strip() for x in line.strip().split(' ') if x.strip() != '']

            if len(line_data) < 6:
                if len(single_structure) > 0:
                    structure_list.append(single_structure)
                    single_structure = []
            else:
                if int(line_data[4]) != 0:
                    base_pair = (int(line_data[0]) - 1, int(line_data[4]) - 1)
                    if base_pair[0] < base_pair[1]:
                        single_structure.append(base_pair)
                
        if len(single_structure) > 0:
            structure_list.append(single_structure)

    sequence_letters = []
    last_index = 0
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            line_data = [x.strip() for x in line.strip().split(' ') if x.strip() != '']

            if len(line_data) < 6:
                continue

            if int(line_data[0]) == last_index + 1:
                sequence_letters.append(line_data[1].upper())
                last_index = int(line_data[0])
            else:
                break

    return structure_list, tuple(sequence_letters)

#converts a structure of basepairs into a structure of helices
def Basepairs_To_Helices(basepair_list):

    if len(basepair_list) == 0:
        return []
    
    basepair_list.sort()

    helix_list = []
    current_helix_start = last_basepair = basepair_list[0]
    current_helix_len = 1
    for basepair in basepair_list[1:]:
        if basepair[0] != last_basepair[0] + 1 or \
                basepair[1] != last_basepair[1] - 1:

            helix_list.append((current_helix_start[0],
                            current_helix_start[1],
                            current_helix_len))
            current_helix_start = last_basepair = basepair
            current_helix_len = 1
        else:
            last_basepair = basepair
            current_helix_len += 1

    helix_list.append((current_helix_start[0],
                    current_helix_start[1],
                    current_helix_len))

    return helix_list

def Helices_To_Basepairs(helix_list, shift=False):
    bps = []
    for (i,j,k) in helix_list:
        if shift:
            bps += [(i+n-1,j-n-1) for n in range(k)]
        else:
            bps += [(i+n,j-n) for n in range(k)]

    return bps

def Helices_To_Helix_Classes(helix_list, sequence):
    result = set()
    for helix in helix_list:
        helix_class = Helix_To_Helix_Class(helix, sequence)

        if helix_class is None:
            print("ERROR: non-canonical basepair present in helix ", helix)
            quit()

        result.add(helix_class)

    return tuple(sorted(list(result)))

def Helices_To_Helix_Class_Dict(helix_list, sequence):
    result = {}
    for helix in helix_list:
        helix_class = Helix_To_Helix_Class(helix, sequence)
        result[helix] = helix_class

    return result

def Helix_Classes_To_Profiles(helix_class_list, featured_helix_classes):
    result = [helix_class for helix_class in helix_class_list 
                if helix_class in featured_helix_classes]
    return tuple(sorted(result))

from Enumerate import get_paired_letters
def is_Valid_Pair(base_a, base_b):
    
    paired_letters = get_paired_letters(base_a)

    if base_b in paired_letters:
        return True
    return False

    '''
    pair = sorted([base_a, base_b])

    if pair[0] == 'A' and pair[1] == 'U':
        return True
    if pair[0] == 'C' and pair[1] == 'G':
        return True
    if pair[0] == 'G' and pair[1] == 'U':
        return True

    return False
    '''
    

def _Helix_To_Helix_Class_impl(helix_tuple, sequence):

    #iterate over basepairs in helix_tuple and verify that they are all valid
    for idx in range(helix_tuple[2]):
        if not is_Valid_Pair(
                sequence[helix_tuple[0] + idx], 
                sequence[helix_tuple[1] - idx]):
            return None

    #iterate from beginning of helix_tuple backwards for as long as possible
    backwards_idx = 0
    while helix_tuple[0] - backwards_idx - 1 >= 0 and\
            helix_tuple[1] + backwards_idx + 1 < len(sequence) and\
            is_Valid_Pair(
                sequence[helix_tuple[0] - backwards_idx - 1],
                sequence[helix_tuple[1] + backwards_idx + 1]):
        backwards_idx += 1

    #iterate from end of helix_tuple forwards for as long as possible
    min_hairpin = 3
    forwards_idx = 0
    while helix_tuple[0] + helix_tuple[2] + forwards_idx < \
            helix_tuple[1] - helix_tuple[2] - forwards_idx - min_hairpin and\
            is_Valid_Pair(
                sequence[helix_tuple[0] + helix_tuple[2] + forwards_idx],
                sequence[helix_tuple[1] - helix_tuple[2] - forwards_idx]):
        forwards_idx += 1

    helix_class = (helix_tuple[0] - backwards_idx + 1,
                    helix_tuple[1] + backwards_idx + 1,
                    helix_tuple[2] + backwards_idx + forwards_idx)

    return helix_class

#helices are 0 indexed, helix classes are 1 indexed in this file
helix_class_dict = {}
def Helix_To_Helix_Class(helix_tuple, sequence):
    
    #check for entry in helix_class_dict[sequence][helix_tuple]
    if sequence in helix_class_dict.keys():
        if helix_tuple in helix_class_dict[sequence].keys():
            return helix_class_dict[sequence][helix_tuple]

    helix_class = _Helix_To_Helix_Class_impl(helix_tuple, sequence)

    #add entry to helix_class_dict
    if not sequence in helix_class_dict.keys():
        helix_class_dict[sequence] = {}
    helix_class_dict[sequence][helix_tuple] = helix_class

    return helix_class

def Read_Featured_Helix_Classes(output_file):
    featured_helix_class_list = []
    featured_helix_class_labels = {}

    with open(output_file, 'r') as f:
        for line in f:
            if not line.startswith("Featured helix"):
                continue

            tokens = line.split()
            
            if not len(tokens) == 9:
                print("Wrong number of tokens in line: {}".format(line))
                continue

            #read (i,j,k) helix
            helix_class = (int(tokens[3]), int(tokens[4]), int(tokens[5]))
            helix_class_label = int(tokens[2].strip(':'))

            featured_helix_class_list.append(helix_class)
            featured_helix_class_labels[helix_class] = helix_class_label

    return set(featured_helix_class_list), featured_helix_class_labels

def Read_All_Helix_Classes(output_file):
    helix_class_list = []
    helix_class_labels = {}

    with open(output_file, 'r') as f:
        for line in f:
            if not line.startswith("Helix"):
                continue

            tokens = line.split()
            
            if not len(tokens) == 12:
                print("Wrong number of tokens in line: {}".format(line))
                continue

            #read (i,j,k) helix
            helix_class = (int(tokens[3]), int(tokens[4]), int(tokens[5]))
            helix_class_label = int(tokens[1])

            helix_class_list.append(helix_class)
            helix_class_labels[helix_class] = helix_class_label

    return set(helix_class_list), helix_class_labels

def stemmable(helix_a, helix_b, gap = (2,2)):

    # ensure gap[0] >= gap[1]
    if gap[0] < gap[1]:
        gap = (gap[1], gap[0])

    assert(gap[1] >= 1)

    i, j, k = helix_a
    m, n, o = helix_b

    y = (j + m - i - n) / 2

    if abs((i + j) - (m + n)) > gap[0]:
        return False

    if -o <= y and y <= k:
        return True

    if -o > y:
        if max(abs(m + o - i), abs(n - o - j)) <= gap[0] and \
           min(abs(m + o - i), abs(n - o - j)) <= gap[1]:
            return True

    if y > k:
        if max(abs(i + k - m), abs(j - k - n)) <= gap[0] and \
           min(abs(i + k - m), abs(j - k - n)) <= gap[1]:
            return True

    return False

def alphabetical_label(x):
    from math import floor

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    digits = []

    if x == 0:
        return letters[0]

    while x > 0:
        digits.append(letters[int(x % len(letters))])
        x = floor(x / len(letters))

    digits.reverse()

    return "".join(digits)

def Build_Edge_List(helix_class_list, gap=(2,2)):
    edge_list = []

    for helix_a in helix_class_list:
        for helix_b in helix_class_list:
            if helix_a == helix_b:
                continue

            if stemmable(helix_a, helix_b, gap):
                edge_list.append((helix_a, helix_b))

    return edge_list

def Find_Stems(featured_helix_class_list, featured_labels, return_diameters=False, gap=None):

    if gap is not None:
        edge_list = Build_Edge_List(featured_helix_class_list, gap)
    else:
        edge_list = Build_Edge_List(featured_helix_class_list)

    G = nx.Graph()
    G.add_edges_from(edge_list)
    G.add_nodes_from(featured_helix_class_list)

    stem_list = list(nx.connected_components(G))
    stem_list.sort(key = lambda l: min(featured_labels[x] for x in l))
    stem_list.sort(key = lambda l: -len(l))

    stem_dict = {}
    diameter_dict = {}
    idx = 0
    for stem in stem_list:
        if len(stem) > 1:
            label = alphabetical_label(idx)
            idx += 1
        else:
            label = str(featured_labels[next(iter(stem))])
        for helix_class in stem:
            stem_dict[helix_class] = label

        if return_diameters:
            diameter_dict[label] = nx.algorithms.distance_measures.diameter(
                G.subgraph(stem))

    if return_diameters:
        return stem_dict, diameter_dict

    return stem_dict

def Find_Stemmable_Class_Dict(helix_class_set):

    res = {}
    for helix_class_a in helix_class_set:
        if helix_class_a not in res:
            res[helix_class_a] = set()

        for helix_class_b in helix_class_set:

            if stemmable(helix_class_a, helix_class_b):
                res[helix_class_a].update([helix_class_b])

    return res
    
def Find_Anchored_Stems(helix_class_list, helix_labels, anchoring_stem_dict):

    reversed_stem_dict = dict()
    local_anchoring_dict = dict()
    for key, value in anchoring_stem_dict.items():
        if value in reversed_stem_dict.keys():
            reversed_stem_dict[value].append(key)
        else:
            reversed_stem_dict[value] = [key]
        local_anchoring_dict[key] = value

    edge_list = Build_Edge_List(helix_class_list)

    G = nx.Graph()
    G.add_edges_from(edge_list)
    G.add_nodes_from(helix_class_list)

    helix_classes_in_stems = set()
    queue = deque()

    #Initialize BFS
    for stem in reversed_stem_dict.values():
        for helix_class in stem:
            queue.append(helix_class)
            helix_classes_in_stems.add(helix_class)

    #Breadth first search
    while len(queue) > 0:
        helix_class = queue.popleft()
        for neighbor in G.neighbors(helix_class):
            if neighbor in helix_classes_in_stems:
                continue

            stem_label = local_anchoring_dict[helix_class]

            reversed_stem_dict[stem_label].append(neighbor)
            queue.append(neighbor)
            helix_classes_in_stems.add(neighbor)
            local_anchoring_dict[neighbor] = stem_label

    #Construct list of all helix_classes currently in stems
    anchored_classes = []
    for value in reversed_stem_dict.values():
        anchored_classes += value

    #Handle connected components without anchoring stem
    remaining_helix_classes = set(helix_class_list) - helix_classes_in_stems
    remaining_edge_list = Build_Edge_List(list(remaining_helix_classes))

    G_remaining = nx.Graph()
    G_remaining.add_edges_from(remaining_edge_list)
    G_remaining.add_nodes_from(remaining_helix_classes)

    remaining_stem_list = list(nx.connected_components(G_remaining))
    remaining_stem_list.sort(key = lambda l: min(helix_labels[helix] for helix in l))
    
    #build output featured stem dictionary
    featured_stem_dict = dict()
    for label, helix_classes in reversed_stem_dict.items():
        for helix_class in helix_classes:
            featured_stem_dict[helix_class] = label

    #Build output stem dictionary
    result_stem_dict = dict()
    for label, helix_classes in reversed_stem_dict.items():
        if len(helix_classes) > 1:# and label.isnumeric():
            label+="*"
        for helix_class in helix_classes:
            result_stem_dict[helix_class] = label

    start_alphabet_index = len(reversed_stem_dict)
    label_idx = start_alphabet_index

    for stem in remaining_stem_list:
        if len(stem) > 1:
            label = alphabetical_label(label_idx)
            label_idx += 1
        else:
            label = str(helix_labels[next(iter(stem))])
        for helix_class in stem:
            result_stem_dict[helix_class] = label

    return result_stem_dict, featured_stem_dict

def Helix_Classes_To_Stems(helix_class_list, stem_dict):
    result = set()
    for helix_class in helix_class_list:
        if helix_class in stem_dict:
            result.add(stem_dict[helix_class])

    result = list(result)
    result.sort()
    return tuple(result)

#Read helix sequence from verbose output of RNAprofile; creates a list of the characters in the sequence
def Read_Sequence(inFile):
    with open(inFile, 'r') as f:
        first_line = f.readline()
        sequence = first_line.split()[4]
        return sequence

def Count_Basepairs(structure):
    #helix is (i,j,k) where k is the length of the helix
    total = sum(helix[2] for helix in structure)
    return total

def Count_Featured_Basepairs(structure, sequence, featured_helix_classes):
    total = 0
    for helix in structure:
        helix_class = Helix_To_Helix_Class(helix, sequence)
        
        if helix_class in featured_helix_classes:
            total += helix[2] #helix is (i,j,k)

    return total

def Count_Stem_Basepairs(structure, sequence, stem_dict):
    return_dict = {}
    for helix in structure:
        helix_class = Helix_To_Helix_Class(helix, sequence)

        if stem_dict[helix_class] not in return_dict:
            return_dict[stem_dict[helix_class]] = 0
        return_dict[stem_dict[helix_class]] += helix[2]

    return return_dict

def Find_Stem_Region(helix_class_list):
    left_start = min(helix_class[0] for helix_class in helix_class_list)
    right_end = max(helix_class[1] for helix_class in helix_class_list)
    left_end = max(helix_class[0] + helix_class[2] - 1 for helix_class in helix_class_list)
    right_start = min(helix_class[1] - helix_class[2] + 1 for helix_class in helix_class_list)

    return (left_start, left_end, right_start, right_end)

import itertools
def Region_To_Basepairs(region):
    left_start, left_end, right_start, right_end = region
    res = list(itertools.product(range(left_start, left_end+1),range(right_start, right_end+1)))
    return res

#helix class list is a list of triplets.
# helix_class_list = [(3, 98, 19), (24, 77, 7), (5, 98, 6)]
def Find_Stem_Width(helix_class_list):
    sum_list = [elem[0] + elem[1] for elem in helix_class_list]
    min_sum = min(sum_list)
    max_sum = max(sum_list)
    return max_sum - min_sum + 1

def Find_Stem_Length(helix_class_list):
    min_dif = min(elem[1] - elem[0] - 2 * elem[2] for elem in helix_class_list)
    max_dif = max(elem[1] - elem[0] for elem in helix_class_list)
    return (max_dif - min_dif) / 2

def Generate_Bracket(helix_structures, helix_class_labels, sequence, method=None):
    import random
    from collections import Counter
    import rna_shape

    bracket_list = []
    for structure in random.sample(helix_structures, min(50, len(helix_structures))):
        local_helix_class_dict = Helices_To_Helix_Class_Dict(structure, sequence)
        local_label_dict = {helix:(helix_class_labels[local_helix_class_dict[helix]] 
                            if local_helix_class_dict[helix] in helix_class_labels else "")
                            for helix in structure}
        
        structure_bracket = rna_shape.find_bracket(structure, local_label_dict)
        bracket_list.append(structure_bracket)

    counts = Counter(bracket_list)

    if method is not None and method == "max":
        most_common = max(bracket_list, key = len)
    else:
        most_common = max(bracket_list, key = counts.get)

    if most_common == "":
        most_common = "[]"

    return most_common

def Cluster_Basepairs(basepairs):
    
    import hdbscan
    import numpy as np

    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(basepairs)
    labels = clusterer.labels_

    result_clusters = [[] for _ in range(np.max(labels) + 2)]

    for basepair, index in zip(basepairs, labels):
        result_clusters[index].append(basepair)

    return result_clusters

def helix_dissimilarity(helix_a, helix_b):
    # ensure y >= 0
    if helix_a[1] + helix_b[0] - helix_a[0] - helix_b[1] < 0:
        helix_a, helix_b = helix_b, helix_a

    i, j, k = helix_a
    m, n, o = helix_b

    y = (j + m - i - n) / 2

    if y < k: #Start of helix b is before end of helix a
        return abs(i - m + j - n)

    return abs(i + k - m) + abs(j - k - n)

def stem_dissimilarity(stem_a, stem_b):
    result = float('inf')
    for helix_class_a in stem_a:
        for helix_class_b in stem_b:
            helix_diss = helix_dissimilarity(helix_class_a, helix_class_b)

            result = min(result, helix_diss)

    return result

def Cluster_Stems(stems):
    import hdbscan
    import numpy as np

    dissimilarity_matrix = np.zeros((len(stems),len(stems)))

    for idx_a, stem_a in enumerate(stems):
        for idx_b, stem_b in enumerate(stems):
            if idx_a == idx_b:
                continue

            dissimilarity_matrix[idx_a, idx_b] = stem_dissimilarity(stem_a, stem_b)

    print(dissimilarity_matrix)

    clusterer = hdbscan.HDBSCAN(min_cluster_size=2,cluster_selection_method='leaf')
    clusterer.fit(dissimilarity_matrix)

    labels = clusterer.labels_

    result_clusters = [[] for _ in range(np.max(labels) + 2)]

    for stem, index in zip(stems, labels):
        result_clusters[index].append(stem)

    return result_clusters

def Cluster_Helix_Classes(helix_classes):
    import hdbscan
    import numpy as np

    dissimilarity_matrix = np.zeros((len(helix_classes),len(helix_classes)))

    for idx_a, helix_a in enumerate(helix_classes):
        for idx_b, helix_b in enumerate(helix_classes):
            if idx_a == idx_b:
                continue

            dissimilarity_matrix[idx_a, idx_b] = helix_dissimilarity(helix_a, helix_b)

    clusterer = hdbscan.HDBSCAN(min_cluster_size=2,cluster_selection_method='leaf')
    clusterer.fit(dissimilarity_matrix)

    labels = clusterer.labels_

    result_clusters = [[] for _ in range(np.max(labels) + 2)]

    for helix, index in zip(helix_classes, labels):
        result_clusters[index].append(helix)

    return result_clusters

def flip_dict(dictionary_in):
    dictionary_out = {}
    for key, value in dictionary_in.items():
        if value in dictionary_out:
            dictionary_out[value].append(key)
        else:
            dictionary_out[value] = [key]

    return dictionary_out

def count_features(structure_list):
    result_counts = {}
    for idx, structure in enumerate(structure_list):
        for feature in structure:
            if feature not in result_counts:
                result_counts[feature] = 0

            result_counts[feature] += 1

    return result_counts

def build_count_dict(structures, basepair_counts, covered_basepair_counts):

    result_dict = {}
    for structure, basepair_count, covered_basepair_count in zip(
            structures, basepair_counts, covered_basepair_counts):
        if structure not in result_dict:
            result_dict[structure] = {}
            result_dict[structure]["num"] = 0
            result_dict[structure]["basepair_count"] = 0
            result_dict[structure]["covered_basepair_count"] = 0

        result_dict[structure]["num"] += 1
        result_dict[structure]["basepair_count"] += basepair_count
        result_dict[structure]["covered_basepair_count"] += covered_basepair_count

    return result_dict

def count_dict_to_sorted_list(count_dict):
    result_list = []

    for key, value in count_dict.items():
        list_elem = value.copy()
        list_elem["structure"] = key
        result_list.append(list_elem)

    result_list.sort(key = lambda elem: -elem["num"])
    return result_list

def Get_File_Prefix(sequenceFile):
    sequence_name = os.path.splitext(os.path.basename(sequenceFile))[0]
    if sequence_name.endswith("_seq"):
        sequence_name = sequence_name[:-4]
    return sequence_name

def Get_Data_Dir(sequence_file, seed):
    data_dir = "./test_files/{}/{}".format(
        Get_File_Prefix(sequence_file),
        seed)
    return data_dir

def Prepare_Sequence_Files(
        sequence_file, 
        RNAStructure_directory, 
        RNAProfiling_directory,
        force_recalculate = False,
        seed = 1):

    data_dir = Get_Data_Dir(sequence_file, seed)
    
    if not os.path.isdir(data_dir) or force_recalculate:
        os.system('bash compute_rna_data.sh {} {} {} {}'.format(
            RNAStructure_directory,
            RNAProfiling_directory,
            sequence_file,
            seed))

def save_results_to_file(filename, sorted_count_list, sequence_name, count_type, total_structures, total_basepairs, clear_file=False):
    if clear_file:
        with open(filename, 'w') as f:
            f.write("sequence, grouping_type, profile_count, structure_prop, basepair_prop, full_basepair_prop\n")
    
    with open(filename, 'a') as f:
        seen_structures = 0
        covered_basepairs = 0
        seen_basepairs = 0

        for idx, row in enumerate(sorted_count_list):
            seen_structures += row["num"]
            covered_basepairs += row["covered_basepair_count"]
            seen_basepairs += row["basepair_count"]

            f.write("{},{},{},{},{},{}\n".format(
                sequence_name,
                count_type,
                idx + 1,
                seen_structures / total_structures,
                covered_basepairs / seen_basepairs,
                covered_basepairs / total_basepairs))

def make_stem_plot(stem_dict, sequence_length, filename = None, featured_stems = None):
    from matplotlib import pyplot as plt
    import draw

    reversed_stem_dict = flip_dict(stem_dict)

    extended_edge_list = []
    extended_node_list = []
    for label, stem in reversed_stem_dict.items():
        #if label in printed_extended_stems:
        extended_edge_list += Build_Edge_List(stem)
        extended_node_list += stem

    extended_helix_edge_graph = nx.Graph()
    extended_helix_edge_graph.add_edges_from(extended_edge_list)
    extended_helix_edge_graph.add_nodes_from(extended_node_list)

    extended_component_labels = {}
    for component in nx.connected_components(extended_helix_edge_graph):
        first_helix = next(iter(component))
        label = stem_dict[first_helix]
        extended_component_labels[tuple(sorted(component))] = label

    fig, ax = plt.subplots()
    draw.plot_helix_class_diagram(
        ax,
        extended_helix_edge_graph,
        featured_stems,
        sequence_length = sequence_length,
        component_labels = extended_component_labels,
        filename=filename)

def group_structures_by_structures(low_level_structures, high_level_structures):
    result_dict = {}
    for low, high in zip(low_level_structures, high_level_structures):
        if high not in result_dict:
            result_dict[high] = []
        result_dict[high].append(low)

    return result_dict

if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) != 4:
        print("Usage: python3 data.py sequence_file.txt path/to/RNAStructure/directory/ path/to/RNAProfiling/executable")
        print("       python3 data.py sequence_folder/ path/to/RNAStructure/directory/ path/to/RNAProfiling/executable")
        exit()

    sequence_files = []
    if os.path.isdir(sys.argv[1]):
        for root, dirs, files in os.walk(sys.argv[1]):
            for single_file in files:
                if single_file.endswith(".txt"):
                    sequence_files.append(os.path.join(root,single_file))
    else:
        sequence_files = [sys.argv[1]]

    rnastructure_dir = sys.argv[2]
    rnaprofiling_executable = sys.argv[3]

    seed_list = [1]

    for sequence in sequence_files:
        for seed in seed_list:
            Prepare_Sequence_Files(
                    sequence, 
                    rnastructure_dir, 
                    rnaprofiling_executable, 
                    seed=seed)

    restart_file = True
    outfile = "results.csv"

    for sequence_file in sequence_files:
        sequence_name = Get_File_Prefix(sequence_file)

        print("\n--------------------------------------------------------------------------------------------")
        print("\nStarting {}".format(sequence_name))

        for seed in seed_list:
            print("\nStarting seed {}\n".format(seed))

            data_dir = Get_Data_Dir(sequence_file, seed)

            rnaprofile_output = data_dir + "/output.txt"
            sample_file = data_dir + "/" + sequence_name + ".ct"

            featured_classes, featured_labels = Read_Featured_Helix_Classes(rnaprofile_output)
            helix_classes, helix_labels = Read_All_Helix_Classes(rnaprofile_output)
            stem_dict = Find_Stems(featured_classes, featured_labels)
            extended_stem_dict, feat_ext_stem_dict \
                = Find_Anchored_Stems(helix_classes, helix_labels, stem_dict)

            basepair_structures, sequence_tuple = Read_CT(sample_file)
            sequence = ''.join(sequence_tuple)
            helix_structures = [Basepairs_To_Helices(structure) 
                                for structure in basepair_structures]

            reversed_stem_dict = flip_dict(stem_dict)
            reversed_extended_stem_dict = flip_dict(extended_stem_dict)

            extended_featured_classes = list(feat_ext_stem_dict.keys())

            #######################################################################################

            helix_class_structures = [Helices_To_Helix_Classes(structure, sequence) 
                                        for structure in helix_structures]
            feat_helix_class_structures = [Helix_Classes_To_Profiles(structure, featured_classes) 
                                        for structure in helix_class_structures]
            feat_stem_structures = [Helix_Classes_To_Stems(structure, stem_dict)
                                        for structure in feat_helix_class_structures]
            feat_ext_stem_structures = [Helix_Classes_To_Stems(structure, feat_ext_stem_dict)
                                        for structure in helix_class_structures]
            extended_stem_structures = [Helix_Classes_To_Stems(structure, extended_stem_dict)
                                        for structure in helix_class_structures]

            stem_counts = count_features(feat_stem_structures)
            extended_stem_counts = count_features(extended_stem_structures)
            extended_count_cutoff = 45
            simplified_extended_stem_dict = {key:value for key, value in extended_stem_dict.items()
                                        if extended_stem_counts[value] > extended_count_cutoff}

            print("| Stem | Region | Frequency | Helix Class List |")
            print("| -- | -- | -- | -- |")
            for key, value in reversed_stem_dict.items():
                region = Find_Stem_Region(value)
                label_list = [featured_labels[x] for x in value]
                print("| {:<3} | {:<20} | {:<4} | {} |".format(
                    key, 
                    str(region), 
                    stem_counts[key], 
                    ", ".join(str(x) for x in label_list)))

            print("")

            print("| Extended Stem | Region | Frequency | Helix Class List |")
            print("| -- | -- | -- | -- |")
            for key, value in reversed_extended_stem_dict.items():
                region = Find_Stem_Region(value)
                label_list = [helix_labels[x] for x in value]

                if extended_stem_counts[key] <= extended_count_cutoff:
                    continue

                print("| {:<3} | {:<20} | {:<4} | {} |".format(
                    key,
                    str(region),
                    extended_stem_counts[key],
                    ", ".join(str(x) for x in label_list)))

            print("")

            '''
            for key, value in sorted(list(reversed_extended_stem_dict.items())):
                region = Find_Stem_Region(value)
                label_list = [helix_labels[x] for x in value]
                print("Extended Stem: {:<3} Region: {:<20} Frequency {:<4} Helix Class List: {}".format(
                    key, str(region), extended_stem_counts[key], label_list))

            print("")
            '''

            basepair_counts = [Count_Basepairs(structure) for structure in helix_structures]
            feat_basepair_counts = [Count_Featured_Basepairs(
                                        structure,
                                        sequence,
                                        featured_classes) for structure in helix_structures]
            feat_ext_basepair_counts = [Count_Featured_Basepairs(
                                        structure,
                                        sequence,
                                        extended_featured_classes) for structure in helix_structures]

            overall_num = len(helix_structures)
            overall_count = sum(basepair_counts)
            overall_feat_count = sum(feat_basepair_counts)
            overall_feat_ext_count = sum(feat_ext_basepair_counts)

            helix_class_count_dict = build_count_dict(
                helix_class_structures, basepair_counts, basepair_counts)
            feat_helix_class_count_dict = build_count_dict(
                feat_helix_class_structures, basepair_counts, feat_basepair_counts)
            feat_stem_count_dict = build_count_dict(
                feat_stem_structures, basepair_counts, feat_basepair_counts)
            feat_ext_stem_count_dict = build_count_dict(
                feat_ext_stem_structures, basepair_counts, feat_ext_basepair_counts)
            extended_stem_count_dict = build_count_dict(
                extended_stem_structures, basepair_counts, basepair_counts)

            feat_helix_class_count_list = count_dict_to_sorted_list(
                feat_helix_class_count_dict)
            feat_stem_count_list = count_dict_to_sorted_list(
                feat_stem_count_dict)
            feat_ext_stem_count_list = count_dict_to_sorted_list(
                feat_ext_stem_count_dict)
            extended_stem_count_list = count_dict_to_sorted_list(
                extended_stem_count_dict)

            stem_structure_helix_structure_dict = group_structures_by_structures(
                helix_structures,
                feat_stem_structures)
            ext_structure_helix_structure_dict = group_structures_by_structures(
                helix_structures,
                extended_stem_structures)

            min_profile_frequency = 20
            print("| Stem Signature | Selected Frequency | Extended Frequency | Basepair % | Extended Basepair % | Basepairs / Base |")
            print("| -- | -- | -- | -- | -- | -- |")
            for row in feat_stem_count_list:
                if row["num"] < min_profile_frequency:
                    continue
                if row["structure"] in feat_ext_stem_count_dict:
                    extended_row = feat_ext_stem_count_dict[row["structure"]]
                else:
                    extended_row = {"num":0, "covered_basepair_count":0, "basepair_count":1}

                helix_structure_list = stem_structure_helix_structure_dict[row["structure"]]
                label = Generate_Bracket(helix_structure_list, stem_dict, sequence)
                print("| {} | {} | {} | {:0.2f} | {:0.2f} | {:0.2f} |".format(
                    label,
                    row["num"],
                    extended_row["num"],
                    row["covered_basepair_count"] / row["basepair_count"],
                    extended_row["covered_basepair_count"] / extended_row["basepair_count"],
                    row["basepair_count"] / row["num"] / len(sequence) * 2))
            
            print("")
            print("| Extended Signature | Stem Signature | Frequency | Basepairs / Base |")
            print("| -- | -- | -- | -- |")
            for row in extended_stem_count_list:
                if row["num"] < min_profile_frequency:
                    continue
                
                helix_structure_list = ext_structure_helix_structure_dict[row["structure"]]
                label = Generate_Bracket(helix_structure_list, extended_stem_dict, sequence)
                stem_label = Generate_Bracket(helix_structure_list, stem_dict, sequence)
                print("| {} | {} | {} | {:0.2f} |".format(
                    label,
                    stem_label,
                    row["num"],
                    row["basepair_count"] / row["num"] / len(sequence) * 2))

            diff_counts = {}
            for stem_a, data_a in feat_stem_count_dict.items():
                for stem_b, data_b in feat_stem_count_dict.items():
                    if stem_a == stem_b:
                        continue

                    additional_a = set(stem_a) - set(stem_b)
                    additional_b = set(stem_b) - set(stem_a)

                    if len(additional_a) > 2 or len(additional_b) > 2:
                        continue

                    diff = tuple(sorted([tuple(sorted(additional_a)), 
                                        tuple(sorted(additional_b))]))

                    if diff not in diff_counts:
                        diff_counts[diff] = {}
                        diff_counts[diff]["num"] = 0
                        diff_counts[diff]["unique_num"] = 0

                    diff_counts[diff]["num"] += data_a["num"] * data_b["num"]
                    diff_counts[diff]["unique_num"] += 1

            sorted_diff_counts = sorted(diff_counts.items(), key = lambda x: -x[1]["num"])
        
            print("")
            print("| Diff | Occurance Count | Sqrt Occurance Count | Unique Occurance Count |")
            print("| -- | -- | -- | -- |")
            for diff, count in sorted_diff_counts:
                if count["unique_num"] / 2 < 2:
                    continue
                if (count["num"] / 2) ** 0.5 < 30:
                    continue

                print("| {} : {} | {:d} | {:0.3f} | {:d} |".format(
                    " ".join(diff[0]), 
                    " ".join(diff[1]), 
                    int(count["num"] / 2),
                    (count["num"] / 2) ** 0.5,
                    int(count["unique_num"] / 2)))

            save_results_to_file(
                outfile, 
                feat_helix_class_count_list,
                sequence_name, "features",
                overall_num, overall_count,
                restart_file)
            restart_file = False
            save_results_to_file(
                outfile, 
                feat_stem_count_list,
                sequence_name, "featured_stems",
                overall_num, overall_count)
            save_results_to_file(
                outfile, 
                feat_ext_stem_count_list,
                sequence_name, "featured_extended_stems",
                overall_num, overall_count)
            save_results_to_file(
                outfile, 
                extended_stem_count_list,
                sequence_name, "extended_stems",
                overall_num, overall_count)

            make_stem_plot(stem_dict, len(sequence), data_dir + "/{}_stem_plot.png".format(sequence_name))
            make_stem_plot(simplified_extended_stem_dict, len(sequence), data_dir + "/{}_extended_stem_plot.png".format(sequence_name), featured_classes)
