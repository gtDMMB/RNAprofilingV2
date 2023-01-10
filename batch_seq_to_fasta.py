
import sys
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage:")
    exit()

in_folder = sys.argv[1]
out_folder = sys.argv[2]

name_set = set()
for filename in Path(in_folder).rglob("*_seq.txt"):
    sequence_name = Path(filename).name[:-8]

    if sequence_name in name_set:
        print("Warning: same name occurs multiple times. Skipping {}".format(
            filename))
        continue

    name_set.add(sequence_name)

    with open(filename, "r") as f:
        sequence = next(f).strip().upper()

    result_filename = sequence_name + ".fasta"
    result_path = Path.joinpath(Path(out_folder), result_filename)

    with open(result_path.resolve(), "w") as f:
        f.write("> {}\n".format(sequence_name))
        f.write("{}\n".format(sequence))

