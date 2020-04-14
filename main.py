import argparse
import shlex
import os
import numpy
import sys
import pickle
import Utilities
import SVTool
import shutil

from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE
from shutil import copyfile
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor


parser = argparse.ArgumentParser(description='Gets the SV consensus.')
parser.add_argument('-f', '--sv_folder', help='Folder containing raw outputs from SV callers.', required=True)
parser.add_argument('-s' '--sample', help='Name of the sample.', required=True)

parser.add_argument('-o', '--output', help='Output file.', default="output.vcf")

parser.add_argument('-m', '--min_overlap', help='Minimum number of SVs in the neighbourhood for the SV to be reported (default 3).', type=int, default=3)

parser.add_argument('-t', '--truth', help='File used for training new model.', required=False)

parser.add_argument('-np', '--no_preprocess', help='Flag used for skipping the preprocessing process - all the preprocessed files should be in temp/ folder.', action="store_true", required=False)

args = parser.parse_args()

sv_files = [f for f in listdir(args.sv_folder) if isfile(join(args.sv_folder, f))]

print(sv_files)

# preprocessing of the files
# problems with no svlen?
if not (args.no_preprocess):
    shutil.rmtree("temp");
    os.mkdir("temp");
    print("Preprocessing files...")

    if (args.truth is not None):
        copyfile(args.truth, "temp/truth.vcf")

    sv_tools = Utilities.preprocessFiles(args.sv_folder)

percDiff = 0.1

X_vector = list()
Y_vector = list()

resulting_svs = list()

consensusId = 1

if (args.truth is None): # load model
    filename = 'pretrained.model'
    loaded_model = pickle.load(open(filename, 'rb'))
i = 0
for svtool in sv_tools:
    if (args.truth is not None):
        if(svtool.tool != "truth"):
            continue
    for sv in svtool.sv_list:
        if(sv.used): continue
        #if(sv.chrom != "chr1"):
        #    continue
        candidates = list()
        candidates.append(sv)
        for svtool2 in sv_tools:
            if(svtool.tool == svtool2.tool):
                continue
            for sv2 in svtool2.sv_list:
                if(sv2.used): continue
                if(sv.chrom != sv2.chrom): # speeds the process up
                    continue
                if(sv2.pos > sv.pos+500): # fix later! it should be dependend on ci or % of svlen
                    break
                if(sv.checkOverlap(sv2)):
                   candidates.append(sv2)
                   break
        if(len(candidates) < 3): # if fewer than 3 then no point in checking it out
            continue
        if (args.truth is not None): # learning phase
            candidates.remove(sv)
            X_vector.append(candidates)
            Y_vector.append(sv)
        else:
            freqDict = Utilities.buildFreqDict(candidates)

            # maybe remove all candidates from svtool once consensus was established based on it?
            (majorityFound, firstMajor) = Utilities.findMajority(sv, freqDict, candidates)
            if(majorityFound):
                newSv = SVariant("consensus", None, firstMajor.chrom, firstMajor.pos, "consensus_"+str(consensusId), firstMajor.ref, firstMajor.end, firstMajor.gt, firstMajor.svlen, firstMajor.svtype, -10, 10, -10, 10)
                consensusId += 1
            else:
                #print("1")
                result = loaded_model.predict(preprocess_X([candidates]))
                pos = result[0]
                end = result[1]
                newSv = SVariant("consensus", None, sv.chrom, int(round(pos)), "consensus_"+str(consensusId), sv.ref, int(round(end)), sv.gt, int(round(pos-end)), sv.svtype, -10, 10, -10, 10)
                consensusId += 1
            resulting_svs.append(newSv)
            Utilities.markUsedCandidates(candidates)

if (args.truth is not None): # learning phase
    X_preprocessed_vector = Utilities.preprocess_X(X_vector)
    Y_preprocessed_vector = Utilities.preprocess_Y(Y_vector)

    X_train, X_test, y_train, y_test = train_test_split(X_preprocessed_vector, Y_preprocessed_vector, test_size=0.33, random_state=42, shuffle=True)
    nn = MLPRegressor(hidden_layer_sizes=(49, 14, 7), solver='lbfgs', max_iter=int(1e8), max_fun=30000, random_state=0)

    print("Creating the model...")

    nn.fit(X_train, y_train)


    nn_score = nn.score(X_test, y_test)
    y_pred = nn.predict(X_test)

    print("Score of model: " + str(nn_score))
    error = abs(y_test-y_pred)
    print("Average abs error (testing set of 10%): " + str(numpy.average(error)) + "std: " + str(numpy.std(error)))

    filename = 'pretrained.model'
    pickle.dump(nn, open(filename, 'wb'))
else:
    header = Utilities.generate_header(args.sample_name)
    with open("output.vcf", 'w') as fout:
        fout.write(header)
        for sv in resulting_svs:
            fout.write(sv.printVcfLine())

    cmd = "cat output.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "+ "\"sort -k1,1V -k2,2n\"" + r"}' > output_sorted.vcf"
    Utilities.execute_command(cmd)

    os.replace("output_sorted.vcf", args.no_preprocess)

#numpy.savetxt("foo.csv", numpy.concatenate((X_test, numpy.vstack((y_test,y_pred)).T), axis=1), delimiter=',', comments="")

# all files are preprocessed now in unified form