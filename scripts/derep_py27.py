#! /usr/bin/env python
__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

'''
derep_seqs.py
Python version of USEARCH's derep command
Fixed the time output to H:M:S format
'''

from optparse import OptionParser
from Bio import SeqIO
import subprocess 
import sys
import time
from datetime import timedelta

start_time = time.time()

# Create the option parser
parser = OptionParser()

# input -i --input
parser.add_option("-i", "--input", action="store", type="string", dest="input_fasta", help="The fasta to dereplicate")

# Grab command line options
(options, args) = parser.parse_args()

ERROR = False
# Check if the input directory exists
if options.input_fasta == None:
    print 'Please enter a fasta to dereplicate.\n\tFor more information use the -h option'
    ERROR = True

# Quit on error
if ERROR == True:
    sys.exit()

# Set variables from command line options
input_fasta = options.input_fasta

# A hash for the unique sequences
sequences = {}

# Number of lines in the file 
records_process = subprocess.Popen(['wc','-l',input_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
result, err = records_process.communicate()
if records_process.returncode != 0:
    raise IOError(err)
records = int(result.strip().split()[0]) /2
print 'There are %s sequences in the fasta' % records

# Read in the fasta
duplicate_time = time.time()
print 'Checking for duplicates'
sequence_count = 0
for seq_record in SeqIO.parse(input_fasta, 'fasta'):
    sequence_count = sequence_count + 1
    sequence = str(seq_record.seq)
    if sequence not in sequences:
        sequences[sequence]=seq_record.id+';size=1;'
    else:
        count = int(sequences[sequence].split('=')[1].rstrip(';')) + 1
        formated_string = sequences[sequence].split('=')[0]+'='+str(count)+';'
        sequences[sequence] = formated_string
        percent_complete = "%.0f%%" % (100.0 * sequence_count/records)
        sys.stdout.write('Percent complete:\t'+percent_complete+'\r')
        sys.stdout.flush()
duplicate_time_complete = time.time()-duplicate_time
print 'It took %s seconds to remove duplicates' % int(duplicate_time_complete)

# Write out our dereplicated fasta
print 'Writing the output file: derep_%s' % input_fasta
output_fasta = 'derep_'+input_fasta
with open(output_fasta, 'w+') as out:
    for sequence in sequences:
        out.write('>'+sequences[sequence]+'\n'+sequence+'\n')

elapsed_time = time.time()-start_time
delta =  timedelta(seconds=elapsed_time)
print 'Total time:\t%s' % str(delta)