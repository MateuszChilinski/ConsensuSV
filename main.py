import argparse
import shlex
import os

from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE
from shutil import copyfile

debug = 1

def generate_header(sample_name):
    fin = open("header", "rt")
    data = fin.read()
    data = data.replace('SAMPLENAME', sample_name)
    fin.close()
    return data+"\n"

def execute_command(cmd):
    if(debug):
        print(cmd)
        process = Popen(cmd, shell=True, stdout=PIPE)
    else:
        process = Popen(cmd, shell=True, stdout=DEVNULL, stderr=STDOUT)
    process.communicate()

def reheader_all(dirFrom, dirTo):
    # create temp header with sample name
    copyfile("header", "header_temp")
    fin = open("header_temp", "rt")
    data = fin.read()
    data = data.replace('SAMPLENAME', args.sample_name)
    fin.close()
    fin = open("header_temp", "wt")
    fin.write(data)
    fin.close()

    # reheader all files
    for file in sv_files:
        cmd = r"bcftools reheader -h header_temp -o " + dirTo + file + " " + dirFrom + file
        if(debug):
            print(cmd)

        process = Popen(cmd, shell=True, stdout=PIPE)
        process.communicate()
        exit_code = process.wait()
    os.remove("header_temp")

class SVariant:
    def __init__(self, line):
        self.parse_line(line)
    def parse_line(self, line):
        values = line.split("\t")
        self.chrom = values[0]
        self.pos = values[1]
        info = values[7].split(";")
        self.gt = values[9]

        self.end = info[0]
        self.svlen = info[1]

        if(self.svlen == "."):
            self.svlen = end-pos

        self.svtype = info[2]
        self.cipos = info[3].split(",")
        self.ciend = info[4].split(",")
        
        if(self.cipos == "."):
            self.cipos = "-10,10" # maybe other values? 0s?
        self.cipos1 = cipos[0]
        self.cipos2 = cipos[1]

        if(self.ciend == "."):
            self.ciend = "-10,10" # maybe other values? 0s?
        self.ciend1 = ciend[0]
        self.ciend2 = ciend[1]

        return
    def print(self):
        print(self.svtype + ": " + self.chrom + " " + self.pos + "(" + self.cipos1 +", " + self.cipos2 + ")" + " - " + self.end + "(" + self.cipos1 +", " + self.cipos2 + ")" + " LEN: " + self.svlen + " GT: " + self.gt)

class SVTool:
    def __init__(self, filename):
        self.parse_file(filename)
    def parse_file(self, filename):
        self.sv_list = list()
        with open(filename) as file:
            line = file.readline()
            while line:
                if not(line.startswith('#')):
                    sv = SVariant(line)
                    sv.print()
                    self.sv_list.add(sv)

parser = argparse.ArgumentParser(description='Gets the SV consensus.')
parser.add_argument('sv_folder', metavar='sv_folder',
                   help='folder consisting the vcf files')
parser.add_argument('sample_name', metavar='sample_name',
                   help='name of the sample')

args = parser.parse_args()

sv_files = [f for f in listdir(args.sv_folder) if isfile(join(args.sv_folder, f))]

print(sv_files)

print("Preprocessing files...")
# preprocessing of the files
# problems with no svlen?
os.mkdir("temp");

reheader_all(args.sv_folder, "temp/")

sv_files = [f for f in listdir("temp/") if isfile(join("temp/", f))]

for file in sv_files:
    # awk -F '\t' '{ $4 = ($4 == "\." ? "N" : $4) } 1' OFS='\t' novoBreak.vcf

    cmd = "sed -i '/:ME:/d' temp/" + file
    execute_command(cmd)

    cmd = "sed -i '/0\/0/d' temp/" + file
    execute_command(cmd)

    cmd = "awk -F " + r"'\t'" + " '{ $4 = ($4 == \"\.\" ? \"N\" : $4) } 1' OFS=" + r"'\t' temp/" + file + " > temp/" + file + "_2"
    execute_command(cmd)
    
    # remove MEI if there are any

    # ensures there are no . in ref
    additional_filters = r"SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND"

    cmd = r"bcftools query -H -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr22,chrX,chrY,chrM -i '(QUAL >= 50 || QUAL = " + "\".\"" + r") && ((SVLEN = " + "\".\"" + r") || (SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50))' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;"+additional_filters+r"\tGT\t[%GT]\n' -o temp/"+file+" temp/"+file+"_2"    
    execute_command(cmd)

    #os.replace("temp/"+file+"_2", "temp/"+file)
    os.remove("temp/"+file+"_2")
    header = generate_header(args.sample_name)
    with open("temp/"+file, 'r') as fin:
        data = fin.read().splitlines(True)
    with open("temp/"+file, 'w') as fout:
        fout.write(header)
        fout.writelines(data[1:])
    svtool = SVTool("temp/"+file)

# all files are preprocessed now in unified form