
def get_paired_letters(letter, allow_ambiguous=True):
    if letter == "A":
        return set(["U","T","N","Y","K","W"])
    if letter == "T":
        return set(["A","N","R","M","W"])
    if letter == "U":
        return set(["A","G","N","R","M","S","W"])
    if letter == "G":
        return set(["C","U","N","Y","M","S"])
    if letter == "C":
        return set(["G","N","R","K","S"])

    if allow_ambiguous:
        if letter == "N":
            return set(["A","G","C","U","T","N","R","Y","K","M","S","W"])
        if letter == "R":
            return set(["C","U","T","N","Y","K","M","S","W"])
        if letter == "Y":
            return set(["A","G","N","R","K","M","S","W"])
        if letter == "K":
            return set(["C","A","N","R","Y","M","S","W"])
        if letter == "M":
            return set(["U","T","G","N","R","Y","K","S","W"])
        if letter == "S":
            return set(["G","C","U","N","R","Y","K","M","S"])
        if letter == "W":
            return set(["A","T","U","N","R","Y","K","M","W"])

        # TODO: there are some other ambiguous nucleotides which 
        # should be included if their absence causes issues

    print("# letter not recognized")
    return set()

def Enumerate(sequence, min_k=1, hairpin_length=3):

    class_list = []
    skip_list = [set() for _ in range(len(sequence))]

    for idx, letter in enumerate(sequence):
        paired_letters = get_paired_letters(letter)

        end_idx = len(sequence) - 1

        while (end_idx >= idx + min_k * 2 - 1 + hairpin_length):

            if (end_idx in skip_list[idx]):
                end_idx -= 1
                continue

            length = 0
            paired_temp = set(paired_letters)

            moving_idx = end_idx
            while ((idx + length + hairpin_length < moving_idx) and
                    (sequence[moving_idx] in paired_temp)):
                skip_list[idx + length].add(moving_idx)

                length += 1
                paired_temp = get_paired_letters(sequence[idx + length])
                moving_idx -= 1

            end_pair = moving_idx + length
            true_length = min(length, int((moving_idx + length - idx + 1) / 2))

            if (true_length >= min_k):
                class_list.append((idx, end_pair, true_length))

            end_idx -= 1

    return class_list

def base_pair_energy(a, b):
    if b < a:
        a, b = b, a

    if a == "A":
        if b == "U":
            return -2

    if a == "C":
        if b == "G":
            return -3

    if a == "G":
        if b == "U":
            return -1

    return 0

def helix_class_energy(i, j, k, sequence):
    energy = sum(
        base_pair_energy(sequence[i + idx], sequence[j - idx]) 
        for idx in range(k))

    return energy

def _generate_helix_class_labels_impl(sequence, min_k=1, hairpin_length=3):
    helix_class_list = Enumerate(sequence, min_k, hairpin_length)

    helix_class_list.sort(key=lambda helix: 
        helix_class_energy(*helix, sequence))

    helix_class_list = [(i + 1, j + 1, k) for i, j, k in helix_class_list]
    helix_class_dict = {helix:str(idx+1) for idx, helix in enumerate(helix_class_list)}

    return helix_class_dict

cached_label_dict = {}
def generate_helix_class_labels(sequence, min_k=1, hairpin_length=3):

    if (sequence, min_k, hairpin_length) in cached_label_dict:
        return cached_label_dict[(sequence, min_k, hairpin_length)]

    result_dict = _generate_helix_class_labels_impl(sequence, min_k, hairpin_length)
    cached_label_dict[(sequence, min_k, hairpin_length)] = result_dict

    return result_dict
    
