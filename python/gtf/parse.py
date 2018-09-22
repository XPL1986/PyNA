# Adapted from: https://gist.github.com/slowkow/8101481

import re
import _pickle as pickle
from interlap import InterLap


def gtf(directory):
    """Parses gtf file from given directory into an interval tree 0 index [start:end) and names"""

    genome, name = dict(), dict()

    with open(directory) as fh:
        for line in fh:
            if not line.startswith('#'):
                _extract(_parse(line), genome, name)

    _save(genome)
    _save_name(name)


def _parse(line):
    """Parses a line in gtf into a dictionary"""

    information = dict()
    fields = line.rstrip().split('\t')

    gtf_header = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
    for i, col in enumerate(gtf_header):
        information[col] = _get_value(fields[i])

    info = [field.split() for field in fields[8].split(";") if len(field.split()) == 2]
    attribute = {key: _get_value(value) for key, value in info}

    fields = ['gene_name', 'exon_number', 'transcript_id', 'transcript_type']
    for key in fields:
        _attempt(key, information, attribute)

    return information


def _get_value(value):
    if not value:
        return None
    value = value.strip('"\'')
    r_comma = re.compile(r'\s*,\s*')
    if ',' in value:
        value = re.split(r_comma, value)

    elif value in ['', '.', 'NA']:
        return ''

    return value


def _attempt(key, result, attribute):
    """Populates pipeline_1 name with given key and attribute name"""
    try:
        result[key] = attribute[key]
    except KeyError:
        result[key] = ''


def _extract(information, genome, name):
    """Creates an intervalTree based on the information of the gtf line"""

    seqname = information["seqname"]
    start, end = int(information["start"]), int(information["end"])
    genome.setdefault(seqname, InterLap()).add((start - 1, end, information))

    gene_name = information["gene_name"]
    name.setdefault(gene_name[0], dict()).setdefault(gene_name, []).append(information)


def _save(genome):
    """Saves the dictionary as intervalTree files"""
    for seqname in genome:
        directory = './pkl/seqname/' + seqname + '.pkl'
        with open(directory, 'wb') as f:
            pickle.dump(genome[seqname], f, -1)
        print("Saved to " + directory)


def _save_name(name):
    """Saves the dictionary as intervalTree files"""
    for letter in name:
        directory = './pkl/name/' + letter + '.pkl'
        with open(directory, 'wb') as f:
            pickle.dump(name[letter], f, -1)
        print("Saved to " + directory)
