#taken and slightly modified from
#https://gist.github.com/larsmans/3116927

"""
Three ways of computing the Hellinger distance between two discrete
probability distributions using NumPy and SciPy.
"""

import numpy as np

_INVSQRT2 = 1 / np.sqrt(2)     # 1/sqrt(2) with default precision np.float64

#hellinger3 in the original gist
def hellinger(p, q):
    return np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) * _INVSQRT2

def hellinger_pandas(p, q):
    #return np.sqrt(np.sum(np.sqrt(p).subtract(np.sqrt(q), fill_value=0) ** 2)) * _INVSQRT2

    dict_p = np.sqrt(p).to_dict()
    dict_q = np.sqrt(q).to_dict()

    diffs = [dict_p.get(key,0) - value for key, value in dict_q.items()]
    diffs += [dict_p[key] for key in dict_p.keys() if key not in dict_q.keys()]

    return np.sqrt(np.sum(np.array(diffs) ** 2)) * _INVSQRT2

def hellinger_numpy(p_key, p_values, q_key, q_values):
    if len(p_key) != len(p_values) or len(q_key) != len(q_values):
        print("WARNING: key and value lengths do not match in hellinger distance calculation")

    p_key = [tuple(key) for key in p_key]
    q_key = [tuple(key) for key in q_key]

    dict_p = {key:value for key, value in zip(p_key, np.sqrt(p_values))}
    dict_q = {key:value for key, value in zip(q_key, np.sqrt(q_values))}
    diffs = [dict_p.get(key,0) - value for key, value in dict_q.items()]
    diffs += [dict_p[key] for key in dict_p.keys() if key not in dict_q.keys()]

    return np.sqrt(np.sum(np.array(diffs) ** 2)) * _INVSQRT2

# same computation as hellinger_numpy
def hellinger_full(p_key, p_values, q_key, q_values):

    p_key = [tuple(key) for key in p_key]
    q_key = [tuple(key) for key in q_key]

    key_list = list(set(p_key + q_key))

    dict_p = {key:sqrt(value) for key, value in zip(p_key, p_values)}
    dict_q = {key:sqrt(value) for key, value in zip(q_key, q_values)}

    result = np.sqrt(np.sum(np.array([
                dict_p.get(key,0)-dict_q.get(key,0) for key in key_list
            ]) ** 2)) * _INVSQRT2

    return result
    
