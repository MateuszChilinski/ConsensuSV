import argparse
import shlex
from subprocess import Popen, PIPE

parser = argparse.ArgumentParser(description='Gets the SV consensus.')
parser.add_argument('sv_folder', metavar='sv_folder',
                   help='folder consisting the vcf files')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')

args = parser.parse_args()
print(args.sv_folder)

# preprocessing of the files
# problems with no svlen?
# bcftools query -i '(SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50)' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND\tGT\t[ %GT]\n' Manta.vcf > MantaxD.vcf
# >50bp, <50?kbp
#cmd = ""
#process = Popen(shlex.split(cmd), stdout=PIPE)
#process.communicate()
#exit_code = process.wait()

