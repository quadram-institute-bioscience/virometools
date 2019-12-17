#!/usr/bin/env python
"""
A script that will only print unique sequences
"""

import sys
import pandas
import argparse 

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

    opt_parser = argparse.ArgumentParser(description='Dereplicate a FASTA file')

    # Collect options
    opt_parser.add_argument('-f', '-i', '--fasta',
                            help='Input file in FASTA format',
                            required=True)

    opt_parser.add_argument('-s', '--by-seq',
                            help='Remove duplicates by sequence (default: by name)')

    opt = opt_parser.parse_args()



    try:
        fp =open(opt.fasta, "r")
    except Exception as e:
        eprint("FATAL ERROR:\n Unable to open FASTA file \"{}\". {}".format(opt.fasta, e))
        exit(1)


    seq_tot = 0
    seq_prn = 0
    if (opt.by_seq):
        # By sequence
        sequences = []
        for name, seq,qual in readfq(fp):
            seq_tot+=1
            if not seq in sequences:
                seq_prn+=1
                print('>{}\n{}'.format(name, seq))
                sequences.append(seq)
    else:
        # By name
        names = []
        for name, seq,qual in readfq(fp):
            seq_tot+=1
            if not name in names:
                seq_prn+=1
                print('>{}\n{}'.format(name, seq))
                names.append(name)
    
    eprint('{}/{} sequences printed'.format(seq_prn, seq_tot)) 


