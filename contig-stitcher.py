import argparse
import sys
import os

###Start of main code
##Beginning with the main argument parser

parser = argparse.ArgumentParser()
parser.add_argument('reference', help = "The main .fa/ .fna file to be concatenated")
parser.add_argument('-l', '--length', help = "The max length of each pseudocontig. Default is 1 Mbp (1,000,000)")
parser.add_argument('-v', '--verbose', help = "Print additional commentary", action = 'store_true')
parser.add_argument('-n','--numbers', help = "Use only numbers for pseudocontig labels", action = 'store_true')

try:
    args = parser.parse_args()
except SystemExit:
    print("No input recognised \nType 'python contig-concat.py -h' for help")
    sys.exit()

args = parser.parse_args()

if args.length:                         #This sets either the user-defined or default length
    try:
        int(args.length)
    except ValueError:
        print("ERROR: Cannot parse length '%s'. Please enter a valid integer\nType 'python contig-concat.py -h' for help" % str(args.length))
        sys.exit()
    maxlength = int(args.length)
else:
    maxlength = 1000000

if args.verbose:                #This sets the boolean for printing additional comments
    commenting = True
else:
    commenting = False

if args.numbers:                #This argument lets you just print numbers as the header for the pseudocontigs
    ps_prefix = ""              #This sets the prefix for the pseudocontig labels in the output file
    print("Suppressing pseudocontig labels")
else:
    ps_prefix = "pseudo_contig_"

try:                            #Can the file be opened successfully?
    open(args.reference)
except IOError:
    print("Error: cannot open %s" % (str(args.reference)))
    sys.exit()

inFile = args.reference
file_root = inFile.split('.')[0]     #Prepare a file root for the output pseudocontig and bed files

ngap = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                #ngap is the series of N's that will be used to separate contigs in the pseudocontigs
fline = 80
                #fline specifies how many characters printer per line in the output fasta file
ngap_length = len(ngap)
                #This is just giving us a numerical length for ngap
ps_number = 1
                #This is the number of pseudocontigs being generated. Is incremented later if necessary

##Step 1: Reads in the supplied file, stores all contigs as an unorder dictionary

print("Step 1/3: Reading in data from %s..." % (str(inFile)))

total_length, contig_count = 0, 0   #total_length = The total summed length of all contigs read in
                                    #contig_count = The count of all contigs read in
contigs = {}                        #A dictionary that will store all of the individual contig sequences
name = ""
with open(inFile) as r:
    for line in r:
        new_contig = False
        line = line.strip()
        if ">" in line:
            name = line
            new_contig = True
            contig_count += 1
        if new_contig:
            contigs[name] = ""
        else:
            contigs[name] = contigs[name] + line
            total_length += len(line)
if commenting:
    print("Total number of contigs read in from %s: %s" % (inFile, str(contig_count)))
    print("Total number of basepairs read in: %s" % (str(total_length)))

##Step 2: Concatenates contigs from dictionary into new pseudocontig dictionary and writes bed file

print("Step 2/3: Concatenating contigs & printing %s.concat.bed..." % (file_root))

ps_label = ps_prefix + str(ps_number)   #This string will be used to index each pseudocontig TODO: Implement join
pseudocontigs = {}                      #This dictionary will contain the pseudocontigs
pseudocontigs[ps_label] = ""
current_length, bed_start, bed_end = 0, 1, 1
bed_file = file_root + ".concat.bed"
if os.path.exists(bed_file):
    os.remove(bed_file)
outFile = open(bed_file, 'a')
if total_length < maxlength:
    for header in contigs.keys():
        pseudocontigs[ps_label] += contigs[header]
        pseudocontigs[ps_label] += ngap
        bed_end = bed_start + len(contigs[header]) - 1
        bed_string = "%s\t%s:\t%s-%s\n" % (ps_label, header, str(bed_start), str(bed_end))
        outFile.write(bed_string)
        bed_start = bed_end + ngap_length + 1
    print("One pseudocontig written. Total length: %s basepairs" % ((str(len(pseudocontigs[ps_label])))))
else:
    print("Pseudocontig length of %s exceeded. Writing multiple pseudocontigs..." % str(maxlength))
    for header in contigs.keys():
        if current_length > maxlength:
            ps_number += 1
            ps_label = ps_prefix + str(ps_number)
            pseudocontigs[ps_label] = ""
            current_length = 0
            bed_start = 1
            bed_end = 1
        pseudocontigs[ps_label] += contigs[header]
        current_length += len(contigs[header])
        bed_end = bed_start + current_length - 1
        if (len(pseudocontigs[ps_label]) + ngap_length) < maxlength:
            pseudocontigs[ps_label] += ngap
            current_length += ngap_length
        bed_string = "%s\t%s:\t%s-%s\n" % (ps_label, header, str(bed_start), str(bed_end))
        outFile.write(bed_string)
        bed_start = bed_end + ngap_length + 1
    print("Total number of pseudocontigs written: %s" % (str(ps_number)))
outFile.close()
print("Bed file written")

#Step 3: Prints pseudocontigs from pseudocontig library to file

print("Step 3/3: Writing %s.concat.fa..." % (str(file_root)))

pscontig_file = file_root + ".concat.fa"
if os.path.exists(pscontig_file):
    os.remove(pscontig_file)
outFile = open(pscontig_file, 'a')
for n in range(1, (ps_number+1)):           #Using range here to preserve pseudocontig order in output file
    ps_label = ps_prefix + str(n)
    if commenting:
        print("Writing %s" % ps_label)
    outFile.write(">%s\n" % (str(ps_label)))
    for i in range(0, (len(pseudocontigs[ps_label])+fline), fline):
        substr = pseudocontigs[ps_label][i:(i+fline)]
        outFile.write("%s\n" % (substr))
outFile.close()
print("Done")
