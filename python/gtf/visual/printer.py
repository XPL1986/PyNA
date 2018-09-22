# source: https://stackoverflow.com/questions/5084743/how-to-print-pretty-string-output-in-python

preset = {'standard': [("Seqname", "seqname", 7), ("Feature", "feature", 12), ("Start", "start", 10),
                       ("End", "end", 10), ("Strand", "strand", 6), ("Frame", "frame", 5),
                       ("Gene Name", "gene_name", 14), ("Exon Number", "exon_number", 11),
                       ("Transcript ID", "transcript_id", 20), ('Transcript Type', 'transcript_type', 15)]}


class TablePrinter(object):
    """Print a list of dicts as a table"""

    def __init__(self, fmt, sep=' ', ul=None):
        """        
        @param fmt: list of tuple(heading, key, width)
                        heading: str, column label
                        key: dictionary key to value to print
                        width: int, column width in chars
        @param sep: string, separation between columns
        @param ul: string, character to underline column label, or None for no underlining
        """
        fmt = preset[fmt]
        super(TablePrinter, self).__init__()
        self.fmt = str(sep).join('{lb}{0}:{1}{rb}'.format(key, width, lb='{', rb='}') for heading, key, width in fmt)
        self.head = {key: heading for heading, key, width in fmt}
        self.ul = {key: str(ul) * width for heading, key, width in fmt} if ul else None
        self.width = {key: width for heading, key, width in fmt}

    def row(self, data):
        return self.fmt.format(**{k: str(data.get(k, ''))[:w] for k, w in self.width.items()})

    def __call__(self, data_list):
        _r = self.row
        res = [_r(data) for data in data_list]
        res.insert(0, _r(self.head))
        if self.ul:
            res.insert(1, _r(self.ul))
        return '\n'.join(res)
