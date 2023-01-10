
import sys
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python3 bpseq_to_fasta.py in_folder/ out_folder/")
    exit()

inFolder = sys.argv[1]
outFolder = sys.argv[2]

Path(outFolder).mkdir(parents=True, exist_ok=True)

file_list = list(Path(inFolder).rglob("*.bpseq"))

for file_name in file_list:
    sequence_name = Path(file_name).stem
    fasta_file = Path(outFolder).joinpath(sequence_name + ".fasta")

    sequence = []
    with open(file_name, "r") as f:
        for line in f:
            sequence.append(line.strip().split()[1])

    with open(fasta_file, "w") as f:
        f.write("> {}\n".format(sequence_name))
        f.write("".join(sequence))


