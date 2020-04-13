import argparse
import shlex
import os
import numpy
import sys

from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE
from shutil import copyfile
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor

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
    def __init__(self, line, tool):
        self.tool = tool
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
                    sv = SVariant(line, self.tool)
                    if(abs(sv.ciend2-sv.ciend1) > self.max_conf or abs(sv.cipos2-sv.cipos1) > self.max_conf):
                        continue
                    #print(self.tool + " | ", end = '')
                    #sv.print_sv()
                    self.sv_list.append(sv)

def preprocessFiles(folder):
    reheader_all(folder, "temp/")

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
        
        cmd = "cat temp/" + file + r"_2 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "+ "\"sort -k1,1V -k2,2n\"" + r"}' > temp/" + file
        execute_command(cmd)
        
        # remove MEI if there are any

        # ensures there are no . in ref
        additional_filters = r"SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND"

        cmd = r"bcftools query -H -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr22,chrX,chrY,chrM -i '(QUAL >= 50 || QUAL = " + "\".\"" + r") && ((SVLEN = " + "\".\"" + r") || (SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50))' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;"+additional_filters+r"\tGT\t[%GT]\n' -o temp/"+file+"_2 temp/"+file    
        execute_command(cmd)

        #os.replace("temp/"+file+"_2", "temp/"+file)
        os.replace("temp/"+file+"_2", "temp/"+file)
        header = generate_header(args.sample_name)
        with open("temp/"+file, 'r') as fin:
            data = fin.read().splitlines(True)
        with open("temp/"+file, 'w') as fout:
            fout.write(header)
            fout.writelines(data[1:])
        svtool = SVTool("temp/"+file)
        sv_tools.append(svtool)
    return sv_tools

def buildFreqDict(candidates):
    freqDict = dict()
    for candidate in candidates:
        key = str(candidate.pos)+"-"+str(candidate.end)
        if key not in freqDict:
            freqDict[key] = 1
        else:
            freqDict[key] += 1
    return freqDict

def findMajority(sv, freqDict, candidates):
    majorityFound = False
    for key in freqDict:
        if(freqDict[key]/len(candidates) >= 0.7):
            print(sv.chrom + " " + sv.svtype + " " + str(sv.pos) + " - " + str(sv.end))
            majorityFound = True
            break
    return majorityFound

def createSVTable():
    sv_files = [f for f in listdir("temp/") if isfile(join("temp/", f))]

    sv_tools = list()

    for file in sv_files:
        toolname = file.split(".")[0]
        if(toolname == "truth"):
            continue
        sv_tools.append(toolname)
    sv_tools.sort()
    return sv_tools

def preprocess_Y(Y_vector):
    Y_prepr = list()
    for sv in Y_vector:
        Y_prepr.append(sv.pos)
        Y_prepr.append(sv.end)
    return Y_prepr
def preprocess_X(X_vector):
    X_prepr = list()
    sv_all_tools = createSVTable()
    for candidates in X_vector:
        candidatesY_pos = list()
        candidatesY_end = list()
        for tool in sv_all_tools:
            found = False
            for sv in candidates:
                if(tool == sv.tool):
                    candidatesY_pos.append(sv.pos)
                    candidatesY_end.append(sv.end)
                    found = True
                    break
            if(found):
                continue
            candidatesY_pos.append(0) # tool not present
            candidatesY_end.append(0)
        X_prepr.append(candidatesY_pos)
        X_prepr.append(candidatesY_end)
    return X_prepr

parser = argparse.ArgumentParser(description='Gets the SV consensus.')
parser.add_argument('sv_folder', metavar='sv_folder',
                   help='folder consisting the vcf files')
parser.add_argument('sample_name', metavar='sample_name',
                   help='name of the sample')
parser.add_argument('--truth', help='used for training new model', required=False)

args = parser.parse_args()

sv_files = [f for f in listdir(args.sv_folder) if isfile(join(args.sv_folder, f))]

print(sv_files)

print("Preprocessing files...")
# preprocessing of the files
# problems with no svlen?
os.mkdir("temp");

if (args.truth is not None):
    copyfile(args.truth, "temp/truth.vcf")

sv_tools = preprocessFiles(args.sv_folder)

percDiff = 0.1

X_vector = list()
Y_vector = list()

for svtool in sv_tools:
    if (args.truth is not None):
        if(svtool.tool != "truth"):
            continue
    for sv in svtool.sv_list:
        #if(sv.chrom != "chr1"):
        #    continue
        candidates = list()
        candidates.append(sv)
        for svtool2 in sv_tools:
            if(svtool.tool == svtool2.tool):
                continue
            for sv2 in svtool2.sv_list:
                if(sv.chrom != sv2.chrom): # speeds the process up
                    continue
                if(sv2.pos > sv.pos+500): # fix later! it should be dependend on ci or % of svlen
                    break
                if(sv.checkOverlap(sv2)):
                   candidates.append(sv2)
                   break
        if(len(candidates) < 3): # if fewer than 3 then no point in checking it out
            continue

        freqDict = buildFreqDict(candidates)

        # maybe remove all candidates from svtool once consensus was established based on it?
        majorityFound = findMajority(sv, freqDict, candidates)
        if(majorityFound):
            continue
        if (args.truth is not None):
            #print("Train NN")
            candidates.remove(sv)
            X_vector.append(candidates)
            Y_vector.append(sv)
        else:
            print("Job for NN")

X_preprocessed_vector = preprocess_X(X_vector)
#print(numpy.array(X_preprocessed_vector))

Y_preprocessed_vector = preprocess_Y(Y_vector)
#print(numpy.array(Y_preprocessed_vector))

X_train, X_test, y_train, y_test = train_test_split(X_preprocessed_vector, Y_preprocessed_vector, test_size=0.33, random_state=42, shuffle=True)
nn = MLPRegressor(hidden_layer_sizes=(50, 20, 10), solver='lbfgs', max_iter=int(1e8), random_state=0)
nn.fit(X_train, y_train)

y_pred = nn.predict(X_test)

numpy.set_printoptions(threshold=sys.maxsize)
#print(y_test)
#print(y_pred)

#print(y_test-y_pred)

print("Average abs error: " + str(numpy.average(abs(y_test-y_pred))))

numpy.savetxt("foo.csv", numpy.concatenate((X_test, numpy.vstack((y_test,y_pred)).T), axis=1), delimiter=',', comments="")

# all files are preprocessed now in unified form