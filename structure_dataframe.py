
import numpy as np
import itertools
from copy import copy

class StructureDataframe:
    def __init__(self, structure_list):
        self._columns = sorted(set(itertools.chain.from_iterable(structure_list)))
        self._column_idx_dict = {feature:idx for idx, feature in enumerate(self._columns)}

        structure_array = np.zeros((len(structure_list), len(self._columns)), dtype=bool)
        for structure_idx, structure in enumerate(structure_list):
            for feature in structure:
                structure_array[structure_idx, self._column_idx_dict[feature]] = True

        self._structures, inverse, self._counts = np.unique(
            structure_array, 
            axis=0, 
            return_counts=True, 
            return_inverse=True)

        count_sort = np.argsort(-self._counts)
        self._structures = self._structures[count_sort]
        self._counts = self._counts[count_sort]

        tmp_index_list = [[] for _ in self._structures]
        for idx, structure_idx in enumerate(inverse):
            tmp_index_list[structure_idx].append(idx)

        self._index_list = [[] for _ in self._structures]
        for idx, old_idx in enumerate(count_sort):
            self._index_list[idx] = tmp_index_list[old_idx]

        self._inverse = inverse
        for idx, indices in enumerate(self._index_list):
            self._inverse[indices] = idx

        self._index_mask = np.ones(self._structures.shape[0], dtype=bool)
        self._column_mask = np.ones(self._structures.shape[1], dtype=bool)

        self._array = None
        self._index = None

    def __iter__(self):
        return iter(self.array)

    @property
    def index(self):
        if self._index is None:
            self._index = list(itertools.chain.from_iterable(
                itertools.compress(self._index_list, self._index_mask)))
        return self._index

    @property
    def counts(self):
        return self._counts[self._index_mask]

    @property
    def columns(self):
        return list(itertools.compress(self._columns, self._column_mask))

    @property
    def column_idx_dict(self):
        return self._column_idx_dict

    def __array__(self):
        if self._array is None:
            self._array = self._structures[self._index_mask][:, self._column_mask]
        return self._array

    def __len__(self):
        return self.shape[0]

    def get_original_array(self):
        return self._structures[self._inverse]

    @property
    def array(self):
        return self.__array__()

    @property
    def shape(self):
        if self._array is None:
            return (np.count_nonzero(self._index_mask), np.count_nonzero(self._column_mask))
        return self._array.shape

    def subset(self, implication_list, frequency_cutoff=0):
        number_at_or_above_cutoff = self._counts.shape[0] - np.searchsorted(self._counts[::-1], frequency_cutoff, side="left") 
        
        res = copy(self)
        res._array = None
        res._index = None

        if len(implication_list) == 0:
            return res

        keep_rows = [self._structures[:number_at_or_above_cutoff,self._column_idx_dict[feature]] == value
            for (feature, value) in implication_list]
        keep_rows.append(self._index_mask[:number_at_or_above_cutoff])
        
        res._index_mask = np.zeros(self._index_mask.shape, dtype=bool)
        res._index_mask[:number_at_or_above_cutoff] = np.all(keep_rows, axis=0)

        res._column_mask = copy(self._column_mask)
        for (feature, value) in implication_list:
            res._column_mask[self._column_idx_dict[feature]] = False

        return res

    def __sub__(self, other):
        if not self._comparable(other):
            print("Structures don't match")
            return None #TODO: Make this raise an error

        res = copy(self)
        res._array = None
        res._index = None

        res._index_mask = np.logical_and(self._index_mask, np.logical_not(other._index_mask))

        return res

    def _comparable(self, other):
        if not isinstance(other, StructureDataframe):
            return False

        # checking if they are the same object is faster than the full equality check
        if self._structures is not other._structures:
            return False
        if self._columns is not other._columns:
            return False
        if self._counts is not other._counts:
            return False
        if self._index_list is not other._index_list:
            return False

        return True

    def __eq__(self, other):
        if not self._comparable(other):
            return False

        if np.any(self._index_mask != other._index_mask):
            return False

        return True

    def __hash__(self):
        return hash(tuple(self._index_mask))

    def __lt__(self, other):
        if not self._comparable(other):
            return False

        self_first_true = np.argmax(self._index_mask)
        other_first_true = np.argmax(other._index_mask)

        if self_first_true < other_first_true:
            return True

        if self_first_true == other_first_true:
            if self._index_mask[self_first_true] < other._index_mask[other_first_true]:
                return True

        return False

    def __gt__(self, other):
        if not self._comparable(other):
            return False

        self_first_true = np.argmax(self._index_mask)
        other_first_true = np.argmax(other._index_mask)

        if self_first_true > other_first_true:
            return True

        if self_first_true == other_first_true:
            if self._index_mask[self_first_true] > other._index_mask[other_first_true]:
                return True

        return False

