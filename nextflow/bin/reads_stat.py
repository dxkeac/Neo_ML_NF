#!/usr/bin/env python

#import os
import sys
import subprocess
import re
import shlex

bam = sys.argv[1]
thread = sys.argv[2]

stat = subprocess.check_output(shlex.split('samtools stat -@' + str(thread) + ' ' + bam)).decode('ascii').split('\n')
for line in stat:
    if re.search( r'raw total sequences:', line, re.I):
        total_reads = line.split('\t')[2]
        print('Total reads: {}'.format(total_reads))
    elif re.search( r'reads mapped:', line, re.I):
        mapped_reads = line.split('\t')[2]
        print('Mapped reads: {} ({:.4%})'.format( mapped_reads, int(mapped_reads)/int(total_reads)))
    elif re.search( r'reads properly paired:', line, re.I):
        properly_paired = line.split('\t')[2]
        print('Properly paired reads: {} ({:.4%})'.format( properly_paired, int(properly_paired)/int(total_reads)))
    elif re.search( r'reads duplicated:', line, re.I):
        duplicated_reads = line.split('\t')[2]
        print('Duplicated reads: {} ({:.4%})'.format(duplicated_reads, int(duplicated_reads)/int(total_reads)))
    elif re.search( r'average length:', line, re.I):
        average_length = line.split('\t')[2]
        print('Average length: {}'.format(average_length))
    elif re.search( r'average quality:', line, re.I):
        average_quality = line.split('\t')[2]
        print('Average quality: {}'.format(average_quality))
    elif re.search( r'insert size average:', line, re.I):
        average_insert_size = line.split('\t')[2]
        print('Average insert size: {}'.format(average_insert_size))
    else:
        pass
