import sys, os

# from http://stackoverflow.com/questions/595305/python-path-of-script
pathname = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(os.path.sep.join([pathname, '..', 'analysis_utils']))
sys.path.append(os.path.sep.join([pathname, '..', 'plotting_utils']))

import ParsingUtils
import AnalysisUtils

import numpy as np

def main():
    sims = os.path.sep.join([pathname, '..', 'sims', 'rsem'])

    rdict = {}
    for i in xrange(1, 21):
        path = os.path.sep.join([sims, str(i)])
        print("reading results from {}".format(path))
        tdf = ParsingUtils.readRSEMTruth(os.path.sep.join([path, "truth.tsv"]) , "_true")
        kdf = ParsingUtils.readKallisto(os.path.sep.join([path, "abundance.tsv"]), "_kallisto")
        sdf = ParsingUtils.readSailfish(os.path.sep.join([path, "quant.sf"]), "_sailfish")
        df = tdf.join(kdf, rsuffix="_K").join(sdf, rsuffix="_S")
        rdict[i] = df

    relDiffs = {} 
    for k,v in rdict.iteritems():
        rds = AnalysisUtils.relDiff("TPM_true", "TPM_sailfish", v, verbose=False)
        rdk = AnalysisUtils.relDiff("TPM_true", "TPM_kallisto", v, verbose=False)
        for method, rd in {"sailfish" : rds, "kallisto" : rdk}.iteritems():
            for summaryName, summaryFunc in {"median" : AnalysisUtils.getMedian, "mean" : AnalysisUtils.getMean}.iteritems(): 
                signedKey = "{}_{}_{}".format(method, summaryName, "signed")
                absKey = "{}_{}_{}".format(method, summaryName, "abs")
                if signedKey in relDiffs:
                    relDiffs[signedKey].append(summaryFunc(rd[0]))
                else:
                    relDiffs[signedKey] = [summaryFunc(rd[0])]
                if absKey in relDiffs:
                    relDiffs[absKey].append(summaryFunc(rd[0].abs()))
                else: 
                    relDiffs[absKey] = [summaryFunc(rd[0].abs())]

    for signedness in ["signed", "abs"]:
        for stat in ["median", "mean"]:
            if stat == "median":
                print("mean of medians of {} relative differences is :\n kallisto: {:0.2f}\n sailfish: {:0.2f}\n".format(
                    signedness, np.mean(relDiffs["kallisto_{}_{}".format(stat, signedness)]),
                    np.mean(relDiffs["sailfish_{}_{}".format(stat, signedness)])))
            elif stat == "mean":
                 print("median of means of {} relative differences is :\n kallisto: {:0.2f}\n sailfish: {:0.2f}\n".format(
                    signedness, np.median(relDiffs["kallisto_{}_{}".format(stat, signedness)]),
                    np.median(relDiffs["sailfish_{}_{}".format(stat, signedness)])))


if __name__ == "__main__":
    main()
