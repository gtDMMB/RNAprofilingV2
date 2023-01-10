
from collections import Counter

def gc_content(sequence):
    letter_counts = Counter(sequence)
    result = (letter_counts["G"] + letter_counts["C"]) / len(sequence)
    return result

def get_family_sequences(rfam_seed_file, family, max_length=None):
    
    with open(rfam_seed_file, "r", encoding="latin-1") as f:
        for line in f:
            if not line.startswith("#=GF AC"):
                continue

            line_family = line.strip().split()[-1]

            if line_family == family:
                break

        return_sequences = []
        for line in f:
            if line.isspace():
                continue
            if line.startswith("# STOCKHOLM"):
                break
            if line.startswith("#"):
                continue
            if line.startswith("//"):
                continue

            new_sequence = line.strip().split()[-1].replace("-","").upper()

            if max_length is None:
                return_sequences.append(new_sequence)
            elif len(new_sequence) <= max_length:
                return_sequences.append(new_sequence)

        return list(set(return_sequences))


if __name__ == "__main__":
    
    import sys
    from pathlib import Path

    import data


    if len(sys.argv) != 4:
        print("Usage python3 augment_dataset.py sequence_folder rfam.seed output_folder")
        exit()

    input_folder = sys.argv[1]
    rfam_seed_file = sys.argv[2]
    output_folder = sys.argv[3]

    sequence_files = sorted(data.get_sequences_in_folder(input_folder))

    small_family_count = 0
    no_root_count = 0

    seen_sequences = set()

    for sequence_file in sequence_files:
        sequence, sequence_name = data.Read_FASTA(sequence_file)
        family = sequence_name[:-2]

        family_sequences = get_family_sequences(rfam_seed_file, family)

        if len(family_sequences) < 5:
            #print("Warning: family {} had {} sequences".format(
            #    family,len(family_sequences)))
            small_family_count += 1

        seq_gc = [gc_content(fam_sequence) for fam_sequence in family_sequences]
        root_seq_gc = gc_content(sequence)

        selected_score_list = [root_seq_gc]#[len(sequence)]
        selected_seq_jdx = []

        include_root = True
        if sequence in seen_sequences:
            include_root = False
        seen_sequences.add(sequence)

        seen_root = False
        for idx in range(4):
            
            max_dist = 0
            max_dist_seq_jdx = -1

            for jdx, score in enumerate(seq_gc):
            #for jdx, score in enumerate(len(seq) for seq in family_sequences):
                if family_sequences[jdx] in seen_sequences:
                    continue

                if family_sequences[jdx] == sequence:
                    seen_root = True
                    continue
                
                dist = min((sel_score - score)**2 for sel_score in selected_score_list)

                if dist >= max_dist:
                    max_dist = dist
                    max_dist_seq_jdx = jdx

            if max_dist_seq_jdx == -1:
                break

            selected_seq_jdx.append(max_dist_seq_jdx)
            selected_score_list.append(seq_gc[max_dist_seq_jdx])#len(family_sequences[max_dist_seq_jdx]))

            seen_sequences.add(family_sequences[max_dist_seq_jdx])

        if not seen_root:
            #print("Did not see root sequence for family {}".format(family))
            no_root_count += 1

        if include_root:
            additional_sequences = [sequence] + [family_sequences[jdx] for jdx in selected_seq_jdx]
        else:
            additional_sequences = [family_sequences[jdx] for jdx in selected_seq_jdx]
            print("skipping root of length ",len(sequence))

        #print("Family {} sequences have scores: {}".format(
        #    family, ", ".join(str(x) for x in selected_score_list)))

        if len(family_sequences) == 0:
            print("Family not in rfam:", family)
            family_sequences += [()]

        #print("Min length: {} Max length: {}".format(
        #    min(len(seq) for seq in family_sequences),
        #    max(len(seq) for seq in family_sequences)))

        for idx, add_seq in enumerate(additional_sequences):
            filename = Path(output_folder) / "{}_{}.fasta".format(family,idx)

            with open(filename.resolve(), "w") as f:
                f.write("> {}_{}\n".format(family,idx))
                f.write("{}\n".format(add_seq))


    print("{} families had fewer than 5 sequences".format(small_family_count))
    print("Did not see root for {} families".format(no_root_count))


