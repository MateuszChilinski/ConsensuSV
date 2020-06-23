import re
import argparse
import os

from subprocess import Popen, PIPE
parser = argparse.ArgumentParser(description='Extracts sample.')

parser.add_argument('-s', '--samples', help='Samples name to compare (comma-separated).', required=True)
parser.add_argument('-o', '--outputs', help='Our files (generated by consensuSV, comma-separated).', required=True)

args = parser.parse_args()



our_callers = set(['lumpy', 'novoBreak', 'GenomeStrip', 'VH', 'wham', 'SVelter', 'Manta', 'Pindel', 'Delly']) # + cnvnator, breakseq, breakdancer
samples = args.samples.split(',')
outputs = args.outputs.split(',')

if os.path.exists('charles_pass'):
    os.remove('charles_pass')
if os.path.exists('charles_pass'):
    os.remove('charles_pass')

i = 0
for sample in samples:
    open('charles_pass','w').writelines([ line for line in open('ALL_Illumina_Integrate_20170206.vcf') if 'PASS' in line or '#' in line])
    open('charles_pass2','w').writelines([ line for line in open('charles_pass') if sample in line or '#' in line])

    full_text = ""
    with open('charles_pass2', 'r') as f:
        for line in f:
            if not('#' in line):
                callers = re.findall(r""+sample+r'HG00512:.{0,20},', line)
                infos = line.split('\t')
                callersp = ""
                callers_list = set()
                for caller in callers:
                    caller_parsed = caller.split(',')[0].split(':')[1]
                    if(caller_parsed in callers_list or caller_parsed not in our_callers): continue
                    callersp += caller_parsed + ","
                    callers_list.add(caller_parsed)
                callersp = callersp[0:-1]
                callers_no = callersp.count(',')+1
                if(callers_no >= 3):
                    full_text += '\t'.join((infos[0], infos[1], infos[2], infos[3], infos[4], infos[5], infos[6], "END="+infos[7].split("END=")[1].split(";")[0] + ";ALGORITHMS="+callersp, "GT", "1/1"+"\n"))
            else:
                full_text += line

    with open("charles_pass_final.vcf", 'w') as fout:
        fout.write(full_text)
    output = outputs[i]
    i += 1
    cmd = "bedtools intersect -wa -header -sorted -f 0.8 -r -a charles_pass_final.vcf -b " + output + " -g ../../merging_new_callings/human.hg19.genome > comparison.vcf"
    print(cmd)
    process = Popen(cmd, shell=True, stdout=PIPE)
    process.communicate()
    cmd = "uniq comparison.vcf > uniq_comparison.vcf"
    process = Popen(cmd, shell=True, stdout=PIPE)
    process.communicate()
    cmd = "grep -vc \"#\" charles_pass_final.vcf"
    process = Popen(cmd, shell=True, stdout=PIPE)
    all_charles = str(process.communicate()[0]).split("'")[1].split(r"\n")[0]
    print(all_charles)
    cmd = "grep -vc \"#\" " + output
    process = Popen(cmd, shell=True, stdout=PIPE)
    all_ours = str(process.communicate()[0]).split("'")[1].split(r"\n")[0]
    cmd = "grep -vc \"#\" uniq_comparison.vcf"
    process = Popen(cmd, shell=True, stdout=PIPE)
    all_intersect = str(process.communicate()[0]).split("'")[1].split(r"\n")[0]
    
    print("== SAMPLE: " + sample + " ==")
    print("")
    print("All detected by us: " + all_ours + " All detected by Charles Lee: " + all_charles + " All common between those two sets: " + all_intersect)
    print("We detect " + str(float(all_intersect)/float(all_charles)*100) + "% of Charles Lee SVs using " + str(float(all_intersect)/float(all_ours)*100) + "% of our set")
    print("")
    os.remove('charles_pass')
    os.remove('charles_pass2')
    os.remove('charles_pass_final.vcf')
    os.remove('comparison.vcf')
    os.remove('uniq_comparison.vcf')
