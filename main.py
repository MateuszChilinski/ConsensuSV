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
        self.pos = int(values[1])
        info = values[7].split(";")
        self.gt = values[9]

        self.end = int(info[0].split("=")[1])
        self.svlen = info[1].split("=")[1]

        if(self.svlen == "."):
            self.svlen = self.end-self.pos
        self.svlen = int(self.svlen)

        self.svtype = self.parse_type(info[2].split("=")[1])
        cipos = info[3].split("=")[1]
        ciend = info[4].split("=")[1]

        if(cipos == "."):
            cipos = "-10,10" # maybe other values? 0s?
        if(ciend == "."):
            ciend = "-10,10" # maybe other values? 0s?

        cipos = cipos.split(",")
        ciend = ciend.split(",")

        self.cipos1 = int(cipos[0])
        self.cipos2 = int(cipos[1])

        self.ciend1 = int(ciend[0])
        self.ciend2 = int(ciend[1])
    def parse_type(self, type):
        if "del" in type.casefold():
            return "DEL"
        if "inv" in type.casefold():
            return "INV"
        if "ins" in type.casefold():
            return "INS"
        if "dup" in type.casefold():
            return "DUP"
        return "UNK"
    def print_sv(self):
        print(self.svtype + ": " + self.chrom + " " + str(self.pos) + "(" + str(self.cipos1) +", " + str(self.cipos2) + ")" + " - " + str(self.end) + "(" + str(self.cipos1) +", " + str(self.cipos2) + ")" + " LEN: " + str(self.svlen) + " GT: " + self.gt)
    def checkOverlap(self, sv2):
        if(self.chrom != sv.chrom):
            return False
        # bear in mind that cipos first coord is negative, hence just addition (example cipos=-10,10)
        minPos1 = self.pos+self.cipos1
        maxPos1 = self.pos+self.cipos2
        minPos2 = sv2.pos+sv2.cipos1
        maxPos2 = sv2.pos+sv2.cipos2

        minEnd1 = self.end+self.ciend1
        maxEnd1 = self.end+self.ciend2
        minEnd2 = sv2.end+self.ciend1
        maxEnd2 = sv2.end+self.ciend2
        if(max(minPos1, minPos2) <= min(maxPos1, maxPos2)):
            if(max(minEnd1, minEnd2) <= min(maxEnd1, maxEnd2)):
                return True
        return False

class SVTool:
    max_conf = 200 # max confidence interval length
    def __init__(self, filename):
        self.tool = filename.split("/")[1].split(".")[0]
        self.parse_file(filename)
    def parse_file(self, filename):
        self.sv_list = list()
        with open(filename) as file:
            for line in file:
                if not(line.startswith('#')):
                    sv = SVariant(line)
                    if(abs(sv.ciend2-sv.ciend1) > self.max_conf or abs(sv.cipos2-sv.cipos1) > self.max_conf):
                        continue
                    #print(self.tool + " | ", end = '')
                    #sv.print_sv()
                    self.sv_list.append(sv)

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

sv_tools = list()

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
    sv_tools.append(svtool)

percDiff = 0.1

for svtool in sv_tools:
    for sv in svtool.sv_list:
        candidates = list()
        candidates.append(sv)
        for svtool2 in sv_tools:
            if(svtool.tool == svtool2.tool):
                continue
            for sv2 in svtool2.sv_list:
                if(sv.checkOverlap(sv2)):
                   candidates.append(sv2)
        print(sv.svtype + " " + str(sv.pos) + " - " + str(sv.end))
        for candidate in candidates:
            # create unified one
            print("\t" + candidate.svtype + " " + str(candidate.pos) + " - " + str(candidate.end))
            # maybe remove all candidates from svtool once consensus was established based on it?
# all files are preprocessed now in unified form