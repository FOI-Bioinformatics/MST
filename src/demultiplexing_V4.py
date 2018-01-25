#!/usr/bin/env python

import subprocess
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser._optionals.title = "Options"
parser.add_argument('-m', nargs = 1, help = 'Mapping file. [Required]', required = True)
parser.add_argument('-r1', nargs = 1, help = 'The R1 file (named *_R1_*.fastq.gz). [Required]', required = True)
parser.add_argument('-r2', nargs = 1, help = 'The R2 file (named *_R2_*.fastq.gz). [Required]', required = True)
parser.add_argument('-i1', nargs = 1, help = 'The I1 file (named *_I1_*.fastq.gz). [Required]', required = True)
parser.add_argument('-i2', nargs = 1, help = 'The I2 file (named *_I2_*.fastq.gz). [Required]', required = True)
parser.add_argument('-p', nargs = 1, help = 'Path to the folder where the scripts are stored. [Required]', required = True)
parser.add_argument('-o', nargs = 1, help = 'A name that is given to the output files. [Required]', required = True)
parser.add_argument('--remove_temp', action ='store_true', help = 'If given, the files that are produced during the script will be removed.')
args = parser.parse_args()
arguments = vars(args)

r1 = ''.join(arguments['r1'])
r2 = ''.join(arguments['r2'])
i1 = ''.join(arguments['i1'])
i2 = ''.join(arguments['i2'])
output_name = ''.join(arguments['o'])
scripts_dir = ''.join(arguments['p'])

if arguments['remove_temp']:
    print("********************************************************\n\nAll intermediate files created in the script will be removed when done.\n\n********************************************************")
    remove = True
else:
	remove = False

I1 = i1.split('/')[-1].strip('.gz')
I2 = i2.split('/')[-1].strip('.gz')


print("Demultiplexing started...")

# Make tmp directory and open log file
if not os.path.isdir("deplex_temporary_files_%s" % output_name):
    os.mkdir("deplex_temporary_files_%s" % output_name)
#    log_file = open("deplex_temporary_files_%s/deplex.log" % output_name)
# else:
# 	wtih open("deplex_temporary_files_%s/deplex.log" % output_name) as log_file:
# 		for line in log_file:

#cmd = "mkdir deplex_temporary_files_%s" % output_name
#subprocess.call(cmd, shell = True)

## Unpack and move index files to current folder
child_processes = []
if not os.path.isfile("deplex_temporary_files_%s/%s" % (output_name, I1)):
    cmd1 = "gunzip -c %s > deplex_temporary_files_%s/%s" % (i1, output_name, I1)
    p = subprocess.Popen(cmd1, shell = True)
    child_processes.append(p)
else:
    print("Skipping unzip I1")
if not os.path.isfile("deplex_temporary_files_%s/%s" % (output_name, I1)):
    cmd2 = "gunzip -c %s > deplex_temporary_files_%s/%s" % (i2, output_name, I2)
    p = subprocess.Popen(cmd2, shell = True)
    child_processes.append(p)
else:
    print("Skipping unzip I2")
for cp in child_processes:
    cp.wait()

if not os.path.isdir("depldeplex_temporary_files_%s/check_map/"):
    cmd = "validate_mapping_file.py -m %s -o deplex_temporary_files_%s/check_map" % (''.join(arguments['m']), output_name)
    subprocess.call(cmd, shell = True)
else:
    print("Skipping validate_mapping_file")

mapping_root, ext = os.path.splitext(os.path.split(''.join(arguments['m']))[1])
print(mapping_root)
print("deplex_temporary_files_%s/%s_corrected_1.txt" % (output_name, mapping_root))
if not os.path.isfile("deplex_temporary_files_%s/%s_corrected_1.txt" % (output_name, mapping_root)):
    if not os.path.isfile("deplex_temporary_files_%s/%s_corrected_2.txt" % (output_name, mapping_root)):
        cmd = "%sfix_mappingfile.py deplex_temporary_files_%s/check_map/*_corrected.txt deplex_temporary_files_%s" % (scripts_dir, output_name, output_name)
        subprocess.call(cmd, shell = True)
else:
    print("Skipping fix_mappingfile")
## Demultiplexing and get quality filtering.
child_processes = []
cmd1 = "split_libraries_fastq.py -o deplex_temporary_files_%s/mapping_%s_1 -i %s -b deplex_temporary_files_%s/%s   --rev_comp_mapping_barcodes -m deplex_temporary_files_%s/*_corrected_1.txt  --store_demultiplexed_fastq" % (output_name, output_name, r1, output_name, I1, output_name)
cmd2 = "split_libraries_fastq.py -o deplex_temporary_files_%s/mapping_%s_2 -i %s -b deplex_temporary_files_%s/%s   --rev_comp_mapping_barcodes -m deplex_temporary_files_%s/*_corrected_1.txt  --store_demultiplexed_fastq" % (output_name, output_name, r2, output_name, I1, output_name)
if not os.path.isdir("deplex_temporary_files_%s/mapping_%s_1" % (output_name, output_name)):
    print("In progress: demultiplexing and quality filtering of the R1 data...")
    p = subprocess.Popen(cmd1, shell = True)
    child_processes.append(p)
if not os.path.isdir("deplex_temporary_files_%s/mapping_%s_2" % (output_name, output_name)):
    print("In progress: demultiplexing and quality filtering of the R2 data...")
    p = subprocess.Popen(cmd2, shell = True)
    child_processes.append(p)

for cp in child_processes:
    cp.wait()

## Syncing
cmd = "%ssyncsort_fq deplex_temporary_files_%s/mapping_%s_1/seqs.fastq deplex_temporary_files_%s/mapping_%s_2/seqs.fastq deplex_temporary_files_%s" % (scripts_dir, output_name, output_name, output_name, output_name, output_name)
subprocess.call(cmd, shell = True)

## Cutadapt
## 16S GGATTAGATACCCBDGTAGTC CTGCWGCCNCCCGTAGG Since there are variation among the primers it can be helpful to set the -e (allowed error rate) parameter a little less stringent. default = 0.1
child_processes = []
cmd = "cutadapt -a ATTAGATACCCTAGTAGTCC deplex_temporary_files_%s/mapping_%s_1/seqs.fastq_paired.fq -o deplex_temporary_files_%s/%s_only_R1_clean.fq -e 0.2" % (output_name, output_name, output_name, output_name)
p = subprocess.Popen(cmd, shell = True)
child_processes.append(p)
cmd = "cutadapt -a TTACCGCGGCTGCTGTCAC deplex_temporary_files_%s/mapping_%s_2/seqs.fastq_paired.fq -o deplex_temporary_files_%s/%s_only_R2_clean.fq -e 0.2" % (output_name, output_name, output_name, output_name)
p = subprocess.Popen(cmd, shell = True)
child_processes.append(p)

for cp in child_processes:
    p.wait()

## Flash
cmd = "flash -m 10 -M 250 -r 250 -f 253 -t 12 -o deplex_temporary_files_%s/%s deplex_temporary_files_%s/%s_only_R1_clean.fq  deplex_temporary_files_%s/%s_only_R2_clean.fq" % (output_name, output_name, output_name, output_name, output_name, output_name)
subprocess.call(cmd, shell = True)

## Syncing
cmd = "%ssyncsort_fq_readindex deplex_temporary_files_%s/%s.extendedFrags.fastq deplex_temporary_files_%s/%s deplex_temporary_files_%s/%s" % (scripts_dir, output_name, output_name, output_name, I1, output_name, I2)
subprocess.call(cmd, shell = True)

## Splitting
print("In progress: split_libraries_fastq.py...")
cmd = "split_libraries_fastq.py -o deplex_temporary_files_%s/mapping_%s -i deplex_temporary_files_%s/%s.extendedFrags.fastq_synced.fq -b deplex_temporary_files_%s/%s_synced.fq  --rev_comp_mapping_barcodes -m deplex_temporary_files_%s/*_corrected_1.txt  --store_demultiplexed_fastq --phred_offset 33 -p 0.6" % (output_name, output_name, output_name, output_name, output_name, I1, output_name)
subprocess.call(cmd, shell = True)

cmd = "%ssplit_fastq.py deplex_temporary_files_%s/mapping_%s/seqs.fastq deplex_temporary_files_%s/%s_synced.fq" % (scripts_dir, output_name, output_name, output_name, I2)
out = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True).communicate()

barcodes = out[0].split('\n')[4].split(' ')


for barcode in barcodes:
	print('Analysing ' + barcode)
	cmd = "split_libraries_fastq.py -o deplex_temporary_files_%s/mapping_%s -i deplex_temporary_files_%s/mapping_%s/seqs_%s.fastq -b deplex_temporary_files_%s/mapping_%s/seqs_%s_barcode.fastq   --rev_comp_mapping_barcodes -m deplex_temporary_files_%s/*_corrected_2.txt  --store_demultiplexed_fastq --rev_comp_barcode --phred_offset 33 -p 0.6" % (output_name, barcode, output_name, output_name, barcode, output_name, output_name, barcode, output_name)
	subprocess.call(cmd, shell = True)
	cmd = "chmod 755 -R deplex_temporary_files_%s/mapping_%s" % (output_name, barcode)
	subprocess.call(cmd, shell = True)

cmd = "mkdir deplex_temporary_files_%s/mapping_all" % output_name
subprocess.call(cmd, shell = True)
cmd = "mkdir deplex_temporary_files_%s/mapping_all_fastq" % output_name
subprocess.call(cmd, shell = True)

for barcode in barcodes:
	print('Copying ' + barcode)
	cmd = "cp deplex_temporary_files_%s/mapping_%s/seqs.fna deplex_temporary_files_%s/mapping_all/seqs_%s.fna" % (output_name, barcode, output_name, barcode)
	subprocess.call(cmd, shell = True)
	cmd = "cp deplex_temporary_files_%s/mapping_%s/seqs.fastq deplex_temporary_files_%s/mapping_all_fastq/seqs_%s.fastq" % (output_name, barcode, output_name, barcode)
	subprocess.call(cmd, shell = True)

for barcode in barcodes:
	print('Fixing fasta and fastq file for ' + barcode)
	cmd = "%sfix_header.py deplex_temporary_files_%s/check_map/*_corrected.txt deplex_temporary_files_%s/mapping_all/seqs_%s.fna %s deplex_temporary_files_%s" % (scripts_dir, output_name, output_name, barcode, barcode, output_name)
	subprocess.call(cmd, shell = True)
	cmd = "%sfix_header_fastq.py deplex_temporary_files_%s/check_map/*_corrected.txt deplex_temporary_files_%s/mapping_all_fastq/seqs_%s.fastq %s deplex_temporary_files_%s" % (scripts_dir, output_name, output_name, barcode, barcode, output_name)
	subprocess.call(cmd, shell = True)

## Fix final fasta and fastq files
cmd = "cat deplex_temporary_files_%s/corrected_*.fna > deplex_temporary_files_%s/corrected_tmp.fna" % (output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "cat deplex_temporary_files_%s/corrected_*.fastq > deplex_temporary_files_%s/corrected_tmp.fastq" % (output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "%sfix_ID.py deplex_temporary_files_%s/corrected_tmp.fna deplex_temporary_files_%s" % (scripts_dir, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "%sfix_ID_fastq.py deplex_temporary_files_%s/corrected_tmp.fastq deplex_temporary_files_%s" % (scripts_dir, output_name, output_name)
subprocess.call(cmd, shell = True)

## Move the created files that should be kept to new locations.
cmd = "mkdir %s_demultiplexed_data" % output_name
subprocess.call(cmd, shell = True)
cmd = "mkdir %s_demultiplexed_data/Control" % output_name
subprocess.call(cmd, shell = True)
cmd = "cp deplex_temporary_files_%s/%s.hist* %s_demultiplexed_data/Control/." % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "cp deplex_temporary_files_%s/mapping_%s_1/split_library_log.txt %s_demultiplexed_data/Control/split_library_1_log_read1.txt" % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "cp deplex_temporary_files_%s/mapping_%s_2/split_library_log.txt %s_demultiplexed_data/Control/split_library_1_log_read2.txt" % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "mv deplex_temporary_files_%s/check_map/*_corrected.txt %s_demultiplexed_data/." % (output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "mv deplex_temporary_files_%s/corrected_all.fna %s_demultiplexed_data/%s.fna" % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "mv deplex_temporary_files_%s/corrected_all.fastq %s_demultiplexed_data/%s.fastq" % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)


## Remove the files that are created during the script.
if remove == True:
    cmd = "rm -R deplex_temporary_files_%s" % output_name
    subprocess.call(cmd, shell = True)

cmd = "cat %s_demultiplexed_data/*_corrected.txt | cut -f 1" % output_name
complete = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True).communicate()
all_samples = filter(None, complete[0].split('\n')[1:])

f = open('%s_demultiplexed_data/%s_number_of_reads.txt' % (output_name, output_name), 'w')
f.write('%-20s\t\t%10s\n' % ('#SampleID', 'Number of reads'))
for sampleID in all_samples:
    cmd = "grep %s_ %s_demultiplexed_data/%s.fna | wc -l" % (sampleID, output_name, output_name)
    line = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True).communicate()
    f.write('%-20s \t\t %10s\n' % (sampleID, line[0].split('\n')[0]))
f.close()
'''
## Perform chimera checking with usearch61.
cmd = "identify_chimeric_seqs.py -i %s_demultiplexed_data/%s.fna -o %s_demultiplexed_data/usearch_checked_chimeras_gg/ -m usearch61 -r /mnt/powervault/moaham/QIIME/gg_13_8/rep_set/97_otus.fasta" % (output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
cmd = "filter_fasta.py -f %s_demultiplexed_data/%s.fna -o %s_demultiplexed_data/%s_without_chimeras_gg_based.fna -s %s_demultiplexed_data/usearch_checked_chimeras_gg/non_chimeras.txt" % (output_name, output_name, output_name, output_name, output_name)
subprocess.call(cmd, shell = True)
'''
print("********************************************************\n\nDemultiplexing is done, the fasta file is now available!\n\n********************************************************")

