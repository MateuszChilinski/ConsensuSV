import shlex
import os
import numpy
import sys
import pickle
import utilities
import shutil

from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from input import inputHandling
from SVTools import SVariant

args = inputHandling()



# preprocessing of the files
# problems with no svlen?
if not (args.no_preprocess):
    shutil.rmtree("temp");
    os.mkdir("temp");
    print("Preprocessing files...")

    if (args.truth is not None):
        copyfile(args.truth, "temp/truth.vcf")

    sv_tools = utilities.preprocessFiles(args.sv_folder, args.sample)
else:
    sv_tools = utilities.loadTempFiles(args.sample)
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
            freqDict = utilities.buildFreqDict(candidates)

            # maybe remove all candidates from svtool once consensus was established based on it?
            (majorityFound, firstMajor) = utilities.findMajority(sv, freqDict, candidates)
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
            utilities.markUsedCandidates(candidates)

if (args.truth is not None): # learning phase
    X_preprocessed_vector = utilities.preprocess_X(X_vector)
    Y_preprocessed_vector = utilities.preprocess_Y(Y_vector)

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
    header = utilities.generate_header(args.sample_name)
    with open("output.vcf", 'w') as fout:
        fout.write(header)
        for sv in resulting_svs:
            fout.write(sv.printVcfLine())

    cmd = "cat output.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "+ "\"sort -k1,1V -k2,2n\"" + r"}' > output_sorted.vcf"
    utilities.execute_command(cmd)

    os.replace("output_sorted.vcf", args.no_preprocess)

#numpy.savetxt("foo.csv", numpy.concatenate((X_test, numpy.vstack((y_test,y_pred)).T), axis=1), delimiter=',', comments="")

# all files are preprocessed now in unified form