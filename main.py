import argparse
import shlex
import os
from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE
from shutil import copyfile

def generate_header(sample_name):
    fin = open("header", "rt")
    data = fin.read()
    data = data.replace('SAMPLENAME', sample_name)
    fin.close()
    return data


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


debug = 1

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
    additional_filters = ""

    if file == "fusor.vcf":
        additional_filters = r"SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=0,0;CIEND=0,0"
    else:
        additional_filters = r"SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND"
    cmd = r"bcftools query -H -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr22,chrX,chrY,chrM -i '(SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50)' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;"+additional_filters+r"\tGT\t[ %GT]\n' -o temp/"+file+"_2 temp/"+file    
    if(debug):
        print(cmd)

    process = Popen(cmd, shell=True, stdout=PIPE)
    process.communicate()
    exit_code = process.wait()
    os.replace("temp/"+file+"_2", "temp/"+file)
    
    header = generate_header(args.sample_name)
    with open("temp/"+file, 'r') as fin:
        data = fin.read().splitlines(True)
    with open("temp/"+file, 'w') as fout:
        fout.write(header)
        fout.writelines(data[1:])
