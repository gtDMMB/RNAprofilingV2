
import sys
import time
import itertools
from collections import Counter

import data
import Enumerate
import hasse
from tree import *
from structure_dataframe import StructureDataframe

if len(sys.argv) != 3:
    print("Usage: python3 exp_a.py path/to/sequence/directory output_file_suffix")
    exit()

sequence_folder = sys.argv[1]
sequence_file_list = data.get_sequences_in_folder(sequence_folder)

repeat_count=25

data.Init_RNA_Seed()

output_file_suffix = sys.argv[2]

f_fuzzy = open("fuzzy_summary_data_" + output_file_suffix + ".csv", "w")
f_tree = open("tree_summary_data_" + output_file_suffix + ".csv", "w")
f_stem = open("stem_summary_data_" + output_file_suffix + ".csv", "w")
f_stem_hc = open("stem_hc_data_" + output_file_suffix + ".csv", "w")
f_profile = open("profile_data_" + output_file_suffix + ".csv", "w")

f_hc_coverage = open("hc_coverage_data_" + output_file_suffix + ".csv", "w")

f_hc_coverage.write(",".join(["{}"] * 9).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "hc",
    "hc_in_selected",
    "sel_hc",
    "sel_hc_in_selected"
    ) + "\n")

f_fuzzy.write(",".join(["{}"] * 10).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "number_stems",
    "fuzzy_changed_structures",
    "mean_fuzzy_region_jaccard_index",
    "max_fuzzy_region_jaccard_index",
    "number_overlapping_fuzzy_region_pairs"
    ) + "\n")

#feature_type in helix, helix class, selected helix class, stem, fuzzy stem
f_tree.write(",".join(["{}"] * 24).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "basepairs",
    "feature_type",
    "feature_count",
    "profile_count",
    "profile_basepair_coverage",
    "empty_profile_structure_coverage",
    "selected_profile_count",
    "selected_profile_cutoff_frequency",
    "selected_profile_structure_coverage",
    "selected_profile_basepair_coverage",
    "tree_structure_coverage",
    "tree_edges",
    "hasse_edges",
    "tree_nodes",
    "hasse_nodes",
    "tree_leaves",
    "tree_time",
    "hasse_time",
    "preprocessing_time") + "\n")

f_stem.write(",".join(["{}"]*16).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "stem_label",
    "stem_width",
    "stem_length",
    "stem_diameter",
    "stem_hc_count",
    "stem_frequency",
    "fuzzy_stem_frequency",
    "stem_region_i",
    "stem_region_j",
    "stem_region_k",
    "stem_region_l") + "\n")

f_stem_hc.write(",".join(["{}"]*9).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "stem_label",
    "hc_label",
    "hc_frequency",
    "hc_length"
    ) + "\n")

#for now feature_type = hc
f_profile.write(",".join(["{}"]*9).format(
    "sequence_name",
    "sequence_len",
    "sequence_hash",
    "repeat_idx",
    "structures",
    "feature_type",
    "profile_idx",
    "profile_frequency",
    "feature_label") + "\n")

structure_count = 1000

for sequence_file, repeat_idx in itertools.product(sequence_file_list, list(range(repeat_count))):
    
    print("# " + str(sequence_file) + " --- " + str(repeat_idx))

    data_dict = data.load_sample_sequence(
        sequence_file, 
        repeat_idxs = [repeat_idx],
        structure_count=structure_count,
        cache_folder="cached_structures")

    sequence_name = data_dict["name"]
    sequence = data_dict["sequence"]
    sequence_hash = data_dict["hash"]

    preprocessing_times = {}

    start_time = time.time()

    helix_labels = Enumerate.generate_helix_class_labels(sequence, min_k=1, hairpin_length=3)
    reversed_label_dict = data.flip_dict(helix_labels)

    basepair_structures = next(iter(data_dict["structures"].values()))[0]

    total_basepairs = sum(len(structure) for structure in basepair_structures)
    all_basepairs = set(itertools.chain.from_iterable(basepair_structures))

    helix_structures = [data.Basepairs_To_Helices(structure) 
                        for structure in basepair_structures]
    helix_class_structures = [data.Helices_To_Helix_Classes(structure, sequence) 
                        for structure in helix_structures]

    featured_classes = get_Featured_Helix_Classes(helix_class_structures)

    featured_labels = {hc:helix_labels[hc] 
            for hc in featured_classes}
    reversed_featured_class_dict = data.flip_dict(featured_labels)

    stem_dict, diameter_dict = data.Find_Stems(featured_classes, featured_labels, return_diameters = True)
    reversed_stem_dict = data.flip_dict(stem_dict)

    feat_helix_class_structures = [data.Helix_Classes_To_Profiles(structure, featured_classes) 
                                for structure in helix_class_structures]
    feat_helix_class_structures_lab = [data.Helix_Classes_To_Stems(structure, featured_labels)
                                for structure in helix_class_structures]

    helix_to_basepair_dict = {helix:data.Helices_To_Basepairs([helix]) for helix in set(itertools.chain.from_iterable(helix_structures))}
    helix_class_to_basepair_dict = {helix:data.Helices_To_Basepairs([(helix)], True) for helix in set(itertools.chain.from_iterable(helix_class_structures))}
    helix_class_lab_to_basepair_dict = {label:data.Helices_To_Basepairs(reversed_featured_class_dict[label], True) for label in reversed_featured_class_dict.keys()}

    preprocessing_times["selected_helix_class"] = time.time() - start_time

    feat_stem_structures = [data.Helix_Classes_To_Stems(structure, stem_dict)
                                for structure in feat_helix_class_structures]

    stem_to_basepair_dict = {label:data.Helices_To_Basepairs(reversed_stem_dict[label], True) for label in reversed_stem_dict.keys()}

    preprocessing_times["stem"] = time.time() - start_time

    fuzzy_stem_structures, fuzzy_stem_to_basepair_dict = Fuzz_Stem_Structures(
        reversed_stem_dict, feat_stem_structures, basepair_structures, True)

    preprocessing_times["fuzzy_stem"] = time.time() - start_time

    fuzz_diff_count = sum(1 for fuzz, strict in zip(fuzzy_stem_structures, feat_stem_structures) if fuzz != strict)

    fuzzy_jaccard_index = overlapping_features(fuzzy_stem_to_basepair_dict,all_basepairs)
    
    mean_jaccard_index, max_jaccard_index = 0, 0
    if len(fuzzy_jaccard_index) > 0:
        mean_jaccard_index = sum(value[2] for value in fuzzy_jaccard_index) / len(fuzzy_jaccard_index)
        max_jaccard_index = max(value[2] for value in fuzzy_jaccard_index)

    f_fuzzy.write(",".join(["{}"] * 10).format(
        sequence_name,
        len(sequence),
        sequence_hash,
        repeat_idx,
        structure_count,
        len(fuzzy_stem_to_basepair_dict),
        fuzz_diff_count,
        mean_jaccard_index,
        max_jaccard_index,
        len(["" for value in fuzzy_jaccard_index if value[2] > 0])
        ) + "\n")

    helix_class_counts = data.count_features(helix_class_structures)
    stem_counts = data.count_features(feat_stem_structures)
    fuzzy_stem_counts = data.count_features(fuzzy_stem_structures)

    for feature_type, structure_list, feature_bp_dict in zip(
            ["helix", "helix_class", "selected_helix_class", "stem", "fuzzy_stem"], 
            [helix_structures, helix_class_structures, feat_helix_class_structures_lab, feat_stem_structures, fuzzy_stem_structures],
            [helix_to_basepair_dict, helix_class_to_basepair_dict, helix_class_lab_to_basepair_dict, stem_to_basepair_dict,fuzzy_stem_to_basepair_dict]):

        dataframe = StructureDataframe(structure_list)

        prof_counts = {tuple(key):value for key, value in zip(dataframe, dataframe.counts)}
        selected_count, min_node_freq = data.cutoff_objects_by_entropy(prof_counts.keys(), prof_counts)

        structure_counts = Counter(tuple(sorted(structure)) for structure in structure_list)
        profiles = [struct for struct in 
            sorted(structure_counts.keys(), key = lambda x: -structure_counts[x])
            if structure_counts[struct] >= min_node_freq]

        selected_structure_list = [structure for structure in structure_list if tuple(sorted(structure)) in profiles]
        selected_bp_structure_list = [bp_structure for bp_structure, structure 
            in zip(basepair_structures, structure_list) if tuple(sorted(structure)) in profiles]

        if feature_type == "selected_helix_class":
            selected_hc_structure_list = [hc_structure for hc_structure, structure
                in zip(helix_class_structures, structure_list) if tuple(sorted(structure)) in profiles]

            hc_count = sum(len(struct) for struct in helix_class_structures)
            sel_hc_count = sum(len(struct) for struct in structure_list)

            hc_sel_count = sum(len(struct) for struct in selected_hc_structure_list)
            hc_sel_hc_count = sum(len(struct) for struct in selected_structure_list)

            f_hc_coverage.write(",".join(["{}"] * 9).format(
                sequence_name,
                len(sequence),
                sequence_hash,
                repeat_idx,
                structure_count,
                hc_count,
                hc_sel_count,
                sel_hc_count,
                hc_sel_hc_count
                ) + "\n")

        coverage = sum(count for count in prof_counts.values() if count >= min_node_freq)
        bp_coverage = get_basepair_coverage(structure_list, basepair_structures, feature_bp_dict) 
        selected_structure_bp_coverage = get_basepair_coverage(
            selected_structure_list, selected_bp_structure_list, feature_bp_dict) 

        empty_structure_count = 0
        if () in structure_counts:
            empty_structure_count = structure_counts[()]

        if feature_type not in ["helix", "helix_class"]:
            if len(sequence) < 1000 or feature_type != "selected_helix_class":

                start_time = time.time()
                tree = build_and_clean_tree(dataframe, min_node_freq)
                tree_build_time = time.time() - start_time

                tree_coverage = get_coverage_count(tree)

                tree_edges, tree_nodes = tree.size(), tree.order()
                tree_leaves = len(list(x for x in tree.nodes() if tree.out_degree(x) == 0))
            else:
                tree_time=""
                tree_edges, tree_nodes, tree_leaves = "", "", ""
                tree_coverage = ""

            if len(sequence) < 1000:
                start_time = time.time()
                hasse_diagram = hasse.build_hasse_diagram(profiles)
                hasse_build_time = time.time() - start_time

                hasse_edges, hasse_nodes = hasse_diagram.size(), hasse_diagram.order()
            else:
                hasse_build_time, hasse_edges, hasse_nodes = "", "", ""

        else:
            tree_build_time, hasse_build_time = "", ""
            tree_edges, tree_nodes, tree_leaves = "", "", ""
            hasse_edges, hasse_nodes = "", ""
            tree_coverage = ""

            preprocessing_times[feature_type] = ""

        f_tree.write(",".join(["{}"] * 24).format(
            sequence_name,
            len(sequence),
            sequence_hash,
            repeat_idx,
            structure_count,
            total_basepairs,
            feature_type,
            dataframe.shape[1],
            dataframe.shape[0],
            bp_coverage,
            empty_structure_count,
            len(profiles),
            min_node_freq,
            coverage,
            selected_structure_bp_coverage,
            tree_coverage,
            tree_edges,
            hasse_edges,
            tree_nodes,
            hasse_nodes,
            tree_leaves,
            tree_build_time,
            hasse_build_time,
            preprocessing_times[feature_type]) + "\n")
        
        if feature_type in ["helix", "helix_class"]:
            continue

        for profile_idx, profile in enumerate(profiles):
            profile_feature_string = "|".join(str(feature) for feature in profile)
            f_profile.write(",".join(["{}"]*9).format(
                sequence_name,
                len(sequence),
                sequence_hash,
                repeat_idx,
                structure_count,
                feature_type,
                profile_idx,
                structure_counts[profile],
                profile_feature_string) + "\n")

    for stem in stem_counts.keys():

        region = data.Find_Stem_Region(reversed_stem_dict[stem])
        width = data.Find_Stem_Width(reversed_stem_dict[stem])
        length = data.Find_Stem_Length(reversed_stem_dict[stem])
        diameter = diameter_dict[stem]

        f_stem.write(",".join(["{}"]*16).format(
            sequence_name,
            len(sequence),
            sequence_hash,
            repeat_idx,
            structure_count,
            stem,
            width,
            length,
            diameter,
            len(reversed_stem_dict[stem]),
            stem_counts[stem],
            fuzzy_stem_counts[stem],
            region[0],
            region[1],
            region[2],
            region[3]) + "\n")

        for hc in reversed_stem_dict[stem]:

            f_stem_hc.write(",".join(["{}"]*9).format(
                sequence_name,
                len(sequence),
                sequence_hash,
                repeat_idx,
                structure_count,
                stem,
                helix_labels[hc],
                helix_class_counts[hc],
                hc[2]
                ) + "\n")

f_tree.close()
f_stem.close()
f_stem_hc.close()
