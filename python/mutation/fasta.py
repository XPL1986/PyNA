def load(filename):
    """
    :param filename: fa file for desired chromosome
    :return: list of string of length 50 nucleotides
    """
    directory = './fa/chroms/'
    with open(directory + filename, 'r') as f:
        data = f.read().split()
    f.close()
    return data


def nucleotides(data, start, end):
    """
    :param data: list of string of length 50 nucleotides
    :param start: coordinate 0 index
    :param end: coordinate 0 index
    :return: sequence of nucleotide at chr:start-end
    """

    def nucleotide(index): return data[(index // 50) + 1][index % 50].upper()
    sequence = "".join([nucleotide(index) for index in range(start, end)])
    return sequence

