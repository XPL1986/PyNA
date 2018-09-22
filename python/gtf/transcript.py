# Adapted from: http://seqanswers.com/forums/showthread.php?t=4914

# import HTSeq Only Supported on Mac
import _pickle as pickle


def create_dict(directory, mode):
    """Creates a dictionary with information from given gtf file, exon or CDS only"""
    assert mode in ['exon', 'CDS']

    gtf_file = HTSeq.GFF_Reader(directory, end_included=True)

    temp, data = dict(), dict()

    for feature in gtf_file:
        if feature.type == mode:
            transcript_id = feature.attr['transcript_id']
            gene = temp.setdefault(feature.attr['gene_name'], dict())
            gene.setdefault(transcript_id, []).append(feature)


    for key, value in temp.items():

        gene = data.setdefault(key, dict())
        for transcript_id in sorted(value):
            transcript_length = 0
            for exon in value[transcript_id]:
                transcript_length += exon.iv.length

            gene.setdefault(transcript_id, dict())["length"] = transcript_length
            gene.setdefault(transcript_id, dict())[mode] = [str(feature.iv) for feature in value[transcript_id]]

    _save(data, mode)


def _save(data, mode):
    """Saves the dictionary pickle files"""
    directory = './pkl/dict/' + mode.lower() + '.pkl'
    with open(directory, 'wb') as f:
        pickle.dump(data, f, -1)
    print("Saved to " + directory)

create_dict('/Users/lcwong/Desktop/PyNA/python/gtf/source/gencode.v28.annotation.gtf', 'CDS')