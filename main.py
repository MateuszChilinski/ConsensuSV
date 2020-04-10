import argparse
import shlex
import os
from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description='Gets the SV consensus.')
parser.add_argument('sv_folder', metavar='sv_folder',
                   help='folder consisting the vcf files')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')

args = parser.parse_args()

sv_files = [f for f in listdir(args.sv_folder) if isfile(join(args.sv_folder, f))]

print(sv_files)

print("Preprocessing files...")
# preprocessing of the files
# problems with no svlen?
os.mkdir("temp");
for file in sv_files:
    cmd = r"bcftools query -i '(SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50)' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND\tGT\t[ %GT]\n' "+args.sv_folder+file+" > temp/"+file
    print(cmd)
    process = Popen(cmd, shell=True, stdout=PIPE)
    process.communicate()
    exit_code = process.wait()

