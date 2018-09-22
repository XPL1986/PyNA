import _pickle as pickle
import re
import json
import math

from python.gtf.visual.printer import TablePrinter


def load(directory):
    """Loads pkl file from directory"""

    with open(directory, 'rb') as f:
        data = pickle.load(f)
    return data


def save(data, directory):
    """Saves data into a txt file"""

    with open(directory, "w") as f:
        f.write(json.dumps(data))
    f.close()
    print("Saved to " + directory)


def feature_lookup(gene, feature):
    """Returns every matching feature related to the gene"""

    query = gene.upper()
    try:
        directory = "/Users/lcwong/Desktop/PyNA/python/gtf/pkl/name/" + query[0] + ".pkl"
        data = [line for line in load(directory)[query] if line["feature"] == feature]
        return data

    except KeyError:
        #print(str(query) + " not found.")
        return []


def name_lookup(gene):
    """Returns every feature related to the gene"""

    query = gene.upper()
    try:
        directory = "python/gtf/pkl/name/" + query[0] + ".pkl"
        features = load(directory)[query]
        return features

    except KeyError:
        print(str(query) + " not found.")
        return []


def location_lookup(location):
    """Location formatted as chr:coord with coord in index 1"""

    try:
        chromosome, coordinate = re.split(':', location)
        features = load("python/gtf/pkl/seqname/" + chromosome + ".pkl")
        matches = [data[2] for data in features.find((int(coordinate) - 1, int(coordinate)))]

        if len(matches) != 0:
            return matches

        else:
            distance, gene = math.inf, None
            for interval in features:
                start, end, data = interval[0], interval[1], interval[2]
                test = min([abs(int(start) - int(coordinate)), abs(int(end) - int(coordinate))])
                if distance > test and data['feature'] == 'gene':
                    distance, gene = test, data
            match = [gene]

            print("\nIntergenic. Nearest gene is " + gene["gene_name"] + " with distance(nucleotides) " + str(
                distance) + ".")
            return match

    except (KeyError, ValueError):
        print("Invalid location. Try again with chr:coord")
        return []


def overlap(overlap):
    """Location formatted as chr:start-end with 1 index"""

    try:
        chromosome, start, end = re.split("[:-]", overlap)
        features = load("python/gtf/pkl/seqname/" + chromosome + ".pkl")
        return [data[2] for data in features.find((int(start) - 1, int(end)))]

    except (KeyError, ValueError):
        print("Invalid overlap. Try again with chr:start-end")
        return []


def query(directory):
    """Each row is one gene name"""

    fh = open(directory, 'r')
    data = []
    for line in fh:
        line = line.rstrip().split()
        if line is not []:
            data += [name_lookup(line[0].upper())]
    return data


def view(data):
    print("\n" + TablePrinter('standard', ul='=')(data) + "\n")
