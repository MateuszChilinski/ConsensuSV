import SVariant

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
                    sv = SVariant(self.tool, line)
                    if(abs(sv.ciend2-sv.ciend1) > self.max_conf or abs(sv.cipos2-sv.cipos1) > self.max_conf):
                        continue
                    #print(self.tool + " | ", end = '')
                    #sv.print_sv()
                    self.sv_list.append(sv)