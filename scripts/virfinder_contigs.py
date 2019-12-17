#!/usr/bin/env python

import sys
import pandas
import argparse
import pdb


def eprint(*args, **kwargs):
    """print to STDERR"""
    print(*args, file=sys.stderr, **kwargs)


def verbose(message):
    if opt.verbose:
        eprint(message)

def debug(message):
    if opt.debug:
        eprint('#{}'.format(message))


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


if __name__ == '__main__':

    opt_parser = argparse.ArgumentParser(description='Denoise Illumina cross-talk from OTU tables')

    opt_parser.add_argument('-f', '--contigs-fasta',
                            help='Contigs/Scaffold file in FASTA format',
                            required=True)

    opt_parser.add_argument('-t', '--virfinder-table',
                            help='Virfinder output (CSV) format',
                            required=True)

    opt_parser.add_argument('-p', '--max-p-value',
                            help='Virfinder maximum p-value',
                            type=float,
                            default=0.05)

    opt_parser.add_argument('-s', '--min-score',
                            help='Virfinder minimum score',
                            type=float,
                            default=0.7)

    opt_parser.add_argument('-v', '--verbose',
                            help='Print extra information',
                            action='store_true')

    opt_parser.add_argument('-d', '--debug',
                            help='Print debug information',
                            action='store_true')



    opt = opt_parser.parse_args()

    try:
        data = pandas.read_csv(opt.virfinder_table, sep=',', header=0)
        filtered = data.loc[data['score'] >= opt.min_score].loc[data['pvalue'] <= opt.max_p_value]
    except Exception as e:
        eprint("FATAL ERROR:\nUnable to open table file \"{}\".\n{}".format(opt.virfinder_table, e))
        exit(1)

    try:
        fp =open(opt.contigs_fasta, "r")
    except Exception as e:
        eprint("FATAL ERROR:\n Unable to open contigs file \"{}\". {}".format(opt.contigs_fasta, e))

    verbose("Virfinder output: {}".format(data.shape[0]))
    verbose("Virfinder selected contigs: {}".format(filtered.shape[0]))
    for name, seq, qual in readfq(fp):
        if name in filtered['name'].values:
            print(">{}\n{}".format(name, seq))