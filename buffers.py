"""
Buffer dictionary
"""

class Buffer:
    def __init__(self, types, C_Na, C_Cl, C_tot_list):
        self.types = types              # list of str
        self.C_Na = C_Na                # float
        self.C_Cl = C_Cl                # float
        self.C_tot_list = C_tot_list    # list of float


def buffer(type):
    # Nested dictionary to store buffer parameters
    buffers = {
        "acetate": {
            "pKa_list": [4.76],
            "zA_list": [0, -1],
            "A": 0.5114,
            "b": 0.16
        },
        "bis-tris": {
            "pKa_list": [6.46],
            "zA_list": [1, 0],
            "A": 0.5114,
            "b": 0.1
        },
        "formate": {
            "pKa_list": [3.75],
            "zA_list": [0, -1],
            "A": 0.5114,
            "b": 0.1
        },
        "succinate": {
            "pKa_list": [4.21, 5.64],
            "zA_list": [0, -1, -2],
            "A": 0.5114,
            "b": 0.1
        },
        "glycine": {
            "pKa_list": [2.34, 9.60],
            "zA_list": [1, 0, -1],
            "A": 0.5114,
            "b": 0.1
        },
        "phosphate": {
            "pKa_list": [2.15, 7.21, 12.33],
            "zA_list": [0, -1, -2, -3],
            "A": 0.5114,
            "b": 0.07
        },
        "carbonate": {
            "pKa_list": [6.37, 10.25],
            "zA_list": [0, -1, -2],
            "A": 0.5114,
            "b": 0.1
        },
        "citrate": {
            "pKa_list": [3.13, 4.76, 6.4],
            "zA_list": [0, -1, -2, -3],
            "A": 0.5114,
            "b": 0.1
        },
        "tris": {
            "pKa_list": [8.06],
            "zA_list": [1, 0],
            "A": 0.5114,
            "b": 0.1
        },
        "arginine": {
            "pKa_list": [2.2, 9.00, 12.5],
            "zA_list": [2, 1, 0, -1],
            "A": 0.5114,
            "b": 0.1
        },
        "histidine": {
            "pKa_list": [1.5, 6.1, 9.3],
            "zA_list": [2, 1, 0, -1],
            "A": 0.5114,
            "b": 0.1
        },
    }

    # Look up parameters for buffer
    params = buffers[type]
    pKa_list = params["pKa_list"]
    zA_list = params["zA_list"]
    A = params["A"]  # [(mol/dm^3)^-0.5] Debye-Huckel parameter at 298 K
    b = params["b"]

    return pKa_list, zA_list, A, b


def get_buffer_info(types):
    num_types = len(types)
    pKa_list = [None] * num_types
    zA_list = [None] * num_types
    A_list = [None] * num_types
    b_list = [None] * num_types

    for i, type in enumerate(types):
        pKa, zA, A, b = buffer(type)
        pKa_list[i] = list(pKa)
        zA_list[i] = list(zA)
        A_list[i] = A
        b_list[i] = b

    return pKa_list, zA_list, A_list, b_list
