
if __name__ == "__main__":
    
    import sys
    from hashlib import blake2b

    if len(sys.argv) != 2:
        print("Usage: python3 file.py folder")
        exit()

    folder = sys.argv[1]

    result_dict = {}

    from pathlib import Path
    from itertools import chain
    for fasta_path in chain(
            Path(folder).rglob("*.fasta"),
            Path(folder).rglob("*.FASTA")):

        with open(fasta_path, "r") as f:
            sequence_name = next(f).strip().replace(">","").split(" ")[-1]
            sequence = next(f).strip()

        family_name = Path(fasta_path).parent.name

        h = blake2b(digest_size=6)
        h.update(sequence.encode('ascii','ignore'))
        sequence_hash = h.hexdigest()

        if sequence_name not in result_dict:
            result_dict[sequence_hash] = (sequence_name, family_name)
        else:
            result_dict[sequence_hash] = (sequence_name, "multiple")

    print("sequence_hash,sequence_name,family_name")
    for key, value in result_dict.items():
        print("{},{},{}".format(key, value[0], value[1]))

