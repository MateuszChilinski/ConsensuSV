import re


open('charles_pass','w').writelines([ line for line in open('ALL_Illumina_Integrate_20170206.vcf') if 'PASS' in line or '#' in line])
open('charles_pass2','w').writelines([ line for line in open('charles_pass') if 'HG00512' in line or '#' in line])

full_text = ""

our_callers = set(['lumpy', 'novoBreak', 'GenomeStrip', 'VH', 'wham', 'SVelter', 'Manta', 'Pindel', 'Delly']) # + cnvnator, breakseq, breakdancer

with open('charles_pass2', 'r') as f:
    for line in f:
        if not('#' in line):
            callers = re.findall(r'HG00512:.{0,20},', line)
            infos = line.split('\t')
            callersp = ""
            callers_list = set()
            for caller in callers:
                caller_parsed = caller.split(',')[0].split(':')[1]
                if(caller_parsed in callers_list or caller_parsed not in our_callers): continue
                callersp += caller_parsed + ","
            callersp = callersp[0:-1]
            callers_no = callersp.count(',')+1
            if(callers_no >= 3):
                full_text += '\t'.join((infos[0], infos[1], infos[2], infos[3], infos[4], infos[5], infos[6], "END="+infos[7].split("END=")[1].split(";")[0] + ";ALGORITHMS="+callersp, "GT", "1/1"+"\n"))
        else:
            full_text += line

with open("charles_pass3.vcf", 'w') as fout:
    fout.write(full_text)