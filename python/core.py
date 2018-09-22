import argparse
import sys

import python.gtf.request as request
import python.gtf.transcript as transcript
import python.gtf.parse as parse


class PyNA(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Welcome to PyNA help menu.',
            usage='''lms <mode> <func>
            The most commonly used commands are:
                parse        parses gtf file into dictionary saved as json files
                query       information based on hg38 annotation file

            © Copyright LMS., Inc., 2018. All rights reserved.''')

        parser.add_argument('mode', help='query')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.mode):
            parser.print_help()
            exit(1)
        getattr(self, args.mode)()

    def parse(self):
        parser = argparse.ArgumentParser(
            description='Parses given gtf file')
        parser.add_argument('input', help="directory of gtf file")
        parser.add_argument('-t', '--type', type=str, default=None, help="genome, exon, cds")
        args = parser.parse_args(sys.argv[2:])

        if args.type is None:
            parse.gtf(args.input)
            transcript.create_dict(args.input, "exon")
            transcript.create_dict(args.input, "CDS")

        elif args.type == 'genome':
            parse.gtf(args.input)

        elif args.type == 'exon':
            transcript.create_dict(args.input, args.type)

        elif args.type == 'cds':
            transcript.create_dict(args.input, args.type.upper())

    def query(self):
        parser = argparse.ArgumentParser(
            description='Search for information based on hg38 annotation file')
        parser.add_argument('type', help="name, location, overlap, file")
        args = parser.parse_args(sys.argv[2:3])

        if not hasattr(self, "query_" + args.type):
            parser.print_help()
            exit(1)
        getattr(self, "query_" + args.type)()

    def query_name(self):
        parser = argparse.ArgumentParser(description='Query based on provided name')

        parser.add_argument('names', nargs='+', type=str, help="gene by given name")
        parser.add_argument('-o', '--output', type=str, default='', help="saves as txt file at given directory")
        parser.add_argument('-f', '--feature', type=str, default=None,
                            help="CDS, exon, five_prime_utr, gene, Selenocysteine, \
                                 start_codon, stop_codon, three_prime_utr, transcript")
        args = parser.parse_args(sys.argv[3:])

        data = []
        for name in args.names:
            if args.feature is not None:
                data += (request.feature_lookup(name, args.feature))
            else:
                data += (request.name_lookup(name))
        request.view(data)

        if args.output is not '':
            request.save(data, args.output)

    def query_location(self):
        parser = argparse.ArgumentParser(description='Query based on provided location chr:coord')

        parser.add_argument('location', nargs='+', type=str, help="genes at provided coordinate")
        parser.add_argument('-o', '--output', type=str, default='', help="saves as txt file at given directory")
        args = parser.parse_args(sys.argv[3:])

        data = []
        for location in args.location:
            data += request.location_lookup(location)
        request.view(data)

        if args.output is not '':
            request.save(data, args.output)

    def query_overlap(self):
        parser = argparse.ArgumentParser(description='Query based on overlaps chr:coord-coord')

        parser.add_argument('overlaps', nargs='+', type=str, help="genes at provided intervals")
        parser.add_argument('-o', '--output', type=str, default='', help="saves as txt file at given directory")
        args = parser.parse_args(sys.argv[3:])

        data = []
        for overlap in args.overlaps:
            data += request.overlap(overlap)
        request.view(data)

        if args.output is not '':
            request.save(data, args.output)

    def query_file(self):
        parser = argparse.ArgumentParser(description='Query based on gene name in file')

        parser.add_argument('inputs', nargs='+', type=str, help="file formatted with one gene per line")
        parser.add_argument('-o', '--output', type=str, default='', help="saves as txt file at given directory")
        args = parser.parse_args(sys.argv[3:])

        data = []
        for input in args.inputs:
            data += request.query(input)
        request.view(data)

        if args.output is not '':
            request.save(data, args.output)

from python.mutation import mutation

mutation.snp("chr10", 87933148, "A")
# snp("chr2", 208248418, "T")
