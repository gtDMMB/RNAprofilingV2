
def Find_Sequence_Structure(raw_sequence, bp_list):
    valid_pairs = set([
        ("A", "U"),
        ("C", "G"),
        ("G", "U")
        ])

    sequence_pairs = []
    for basepair in bp_list:
        letter_pair = tuple(sorted([
            raw_sequence[basepair[0]],
            raw_sequence[basepair[1]]]))

        if letter_pair in valid_pairs:
            sequence_pairs.append(basepair)

    return sequence_pairs

def Find_Unaligned_Sequence(raw_sequence, bp_list):
    
    bp_mapping_dict = {}

    new_idx = 0
    sequence = []
    for idx, char in enumerate(raw_sequence):
        if char.isalpha():
            sequence.append(char)
            bp_mapping_dict[idx] = new_idx

            new_idx +=1
            
    res_bp_list = [(bp_mapping_dict[i], bp_mapping_dict[j]) for i,j in bp_list]
        
    return "".join(sequence), res_bp_list

def To_Dot_Bracket(bp_list, sequence_length):
    characters = ["."] * sequence_length

    skipped_pairs = []
    for bp in bp_list:
        characters[bp[0]] = "("
        characters[bp[1]] = ")"

    return "".join(characters)

def To_RNA(sequence):
    return "".join(
        ["T" if (char in {"u", "U"}) else char
        for char in sequence])

if __name__ == "__main__":
    
    import sys    

    import pathlib

    import data
    
    if len(sys.argv) != 3:
        print("Usage: python3 expand_stanford.py filename.stanford outfolder")
        exit()

    infile = sys.argv[1]
    outfolder = sys.argv[2]

    pathlib.Path(outfolder).mkdir(parents=True, exist_ok=True)

    raw_sequence_dict = {}

    with open(infile,'r') as f:
        
        for line in f:
            if line.isspace() or line.startswith("//"):
                continue

            if line.startswith("#=GC SS_cons"):
                consensus_structure = line[12:].strip()
                continue

            if line.startswith("#"):
                continue
            
            line_data = [elem.strip() for elem in line.strip().split()]

            assert(len(line_data) == 2)

            raw_sequence_dict[line_data[0]] = To_RNA(line_data[1].upper())


    bp_list = data.Dot_to_BP(consensus_structure)

    for name, raw_sequence in raw_sequence_dict.items():

        sequence_struct = Find_Sequence_Structure(raw_sequence, bp_list)
        sequence, sequence_struct = Find_Unaligned_Sequence(raw_sequence, sequence_struct)

        cleaned_name = name.replace("/","-")

        fasta_file = pathlib.Path(
            outfolder,
            "{}.fasta".format(cleaned_name)).resolve()
        struct_file = pathlib.Path(
            outfolder,
            "{}.struct".format(cleaned_name)).resolve()

        with open(fasta_file, "w") as f:
            f.write("> {} \n".format(name))
            f.write("{}\n".format(sequence))

        with open(struct_file, "w") as f:
            f.write("> {} \n".format(name))
            f.write("{}\n".format(To_Dot_Bracket(sequence_struct, len(sequence))))


