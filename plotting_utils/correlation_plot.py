import sys
sys.path.append("utilities")
from utilities import ParsingUtils
from utilities import AnalysisUtils

import pandas as pd
import seaborn as sns
import numpy as np
import scipy as sp
from scipy import stats as stats
import matplotlib.pyplot as plt
import argparse
import snakemake

colors = ["dark grey", "greyish", "greyish", "faded green", "amber"]
sns.set_palette(sns.xkcd_palette(colors))

methodcolor = {"Sailfish":"#d7191c", "eXpress":"#fdae61", "Salmon":"#abd9e9", "SalmonAln":"#2c7bb6",
        "Kallisto":"#8c510a", "Salmon (VB)": "#fdae61", "SalmonAln (VB)": "#000000", 
        "Sailfish (Quasi)": "#00ff00", "STAR Alignment": "#0000ff"}

def setPlotProperties():
    # Will these persist properly?

    # Style
    sns.set_style("ticks")

    # Context (poster = large text, and bold)
    sns.set_context("paper")

    # Color palette
    sns.set_palette("pastel")


def makeRelErrorPlot(truthCol, predCol, df, outBase, method, measure):
    plt.cla()
    plt.clf()

    if measure == "tpm":
        outFile = "{}_tpm_rel_error.pdf".format(outBase)
    elif measure == "num_reads":
        outFile = "{}_num_reads_rel_error.pdf".format(outBase)

    # Get the relative difference numbers
    reV = AnalysisUtils.relError(truthCol, predCol, df)
    _ = plt.hist(np.clip(reV, -20.0, 20.0), range=(-20, 20), bins=100, histtype="stepfilled")

    mdText = "mean rel error = {0:.2f}\nmean |rel error| = {1:.2f}\n# data points = {2}".format(
              reV.mean(), reV.abs().mean(), len(reV))

    plt.figtext(0.15, 0.85, mdText)
    
    ax = plt.axes()
    ax.set_ylabel("frequency")
    ax.set_xlabel("relative difference")

    # Get rid of axis spines
    sns.despine()
    plt.savefig(outFile)

def makeRelDiffPlot(truthCol, predCol, df, outBase, method, measure):
    plt.cla()
    plt.clf()

    if measure == "tpm":
        outFile = "{}_tpm_rel_diff.pdf".format(outBase)
    elif measure == "num_reads":
        outFile = "{}_num_reads_rel_diff.pdf".format(outBase)

    # Get the relative difference numbers
    rdV, nzV = AnalysisUtils.relDiff(truthCol, predCol, df)
    try:
        _ = plt.hist(rdV["relDiff"], bins=100, histtype="stepfilled")

        mdText = "mean rel diff = {0:.2f}\nmean |rel diff| = {1:.2f}".format(
                  rdV["relDiff"].mean(), rdV["relDiff"].abs().mean())

        plt.figtext(0.15, 0.85, mdText)
        
        ax = plt.axes()
        ax.set_ylabel("frequency")
        ax.set_xlabel("relative difference")

        # Get rid of axis spines
        sns.despine()
        plt.savefig(outFile)
    except:
        f = open(outFile, 'w')
        f.close()
        print("Error making relative difference plot {}".format(outFile))

def makeCorrPlot(truthCol, predCol, df, outBase, method, measure, logged=True):
    plt.cla()
    plt.clf()

    measureName = "NO MEASURE NAME"
    if measure == "tpm":
        outFile = "{}_tpm_corr.pdf".format(outBase)
        measureName = "TPM"
    elif measure == "num_reads":
        outFile = "{}_num_reads_corr.pdf".format(outBase)
        measureName = "# fragments"


    rv, pv = stats.spearmanr(df[truthCol], df[predCol])
    corrText = "Spearman r = {0:.2f}".format(rv)

    if (logged):
        minLogVal = -2.5
        smallVal = 1e-2
        ax = plt.axes()
        minVal = min(np.log10(df.loc[df[truthCol] > 0, truthCol].min()),
                     np.log10(df.loc[df[predCol] > 0, predCol].min()))
        maxVal = max(np.log10(df.loc[df[truthCol] > 0, truthCol].max()),
                     np.log10(df.loc[df[predCol] > 0, predCol].max()))

        sns.regplot(np.log10(df[truthCol]), np.log10(df[predCol]), df, 
                fit_reg=False, dropna=True, color=[0.7, 0.7, 0.7, 0.2],
                ax=ax)

        ax.set_xlabel("log(True {})".format(measureName))
        ax.set_ylabel("log({} {})".format(method, measureName))

        ax.set_xlim(minVal-0.5, maxVal+0.5) 
        ax.set_ylim(minVal-0.5, maxVal+0.5) 
        plt.figtext(0.15, 0.85, corrText)
        
    else:
        ax = plt.axes()
 
        minVal = min(df[truthCol].min(), df[predCol].min())
        maxVal = max(df[truthCol].max(), df[predCol].max())

        sns.regplot(truthCol, predCol, df, fit_reg=False, 
                    color=[0.7, 0.7, 0.7, 0.2], ax=ax)
 
        ax.set_xlabel("True {}".format(measureName))
        ax.set_ylabel("{} {}".format(method, measureName))

        ax.set_xlim(minVal, maxVal) 
        ax.set_ylim(minVal, maxVal) 
 
        plt.figtext(0.15, 0.85, corrText)

    # Get rid of axis spines
    sns.despine()

    plt.savefig(outFile)


def spearmanCorrelationPlots(D, measure, methods, simIDs, outdir, plotname, setxlim=False, lines=False):
    font = {
            'weight' : 'bold',
            'size'   : 50}

    import matplotlib
    matplotlib.rc('font', **font)
    #sns.palplot(sns.hls_palette(8, l=.3, s=.8))
    #sns.set_palette(sns.hls_palette(8, l=.3, s=.8))
    flatui = ["#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#ffffbf"]
    sns.set_palette(sns.color_palette(flatui))
    sns.set_style('white')
    sns.set_context('poster')
    plt.clf()
    plt.cla()


    md = {}
    for m in methods:
        md[m] = []
        for i in simIDs:
            pc = D["{}_truth{}".format(measure, str(i))].corr(D["{}_{}{}".format(measure, m, str(i))], 'spearman')
            md[m].append(pc)
            print("Method = {}, Pearson corr = {}".format(m, pc))

    ax = None 
    for mn in methods:
        mv = md[mn]
        weights = np.ones_like(mv)/len(mv)
        #ax = sns.distplot(mv, hist_kws={"histtype": "stepfilled"}, kde=False, rug=True, label=mn, ax=ax, norm_hist=False)
        if lines:
            ax = sns.distplot(mv, hist_kws={"histtype": "step", "linewidth": 4, "alpha": 1}, kde=False, rug=False, label=mn, ax=ax, norm_hist=False, color=methodcolor[mn])
        else:
            ax = sns.distplot(mv, hist_kws={"linewidth": 0.1, "alpha": 1}, kde=False, rug=False, label=mn, ax=ax, norm_hist=False, color=methodcolor[mn])

    sns.despine()
    ax = plt.axes()
    #ax.set_xlabel('Spearman correlation')
    #ax.set_ylabel('frequency')
    plt.legend(loc=2)
    if setxlim: 
        plt.xlim(setxlim[0], setxlim[1])
    plt.xlabel("Spearman correlation")
    plt.ylabel("frequency")
    plt.savefig('{}/{}.pdf'.format(outdir, plotname))

def proportionalityCorrelationPlots(D, measure, methods, simIDs, outdir, plotname):
    sns.set_style('white')
    sns.set_context('poster')
    plt.clf()
    plt.cla()

    md = {}
    for m in methods:
        md[m] = []
        for i in simIDs:
            pc = AnalysisUtils.proportionalityCorrelation(
                    "{}_truth{}".format(measure, str(i)), 
                    "{}_{}{}".format(measure, m, str(i)), D)
            md[m].append(pc)

    ax = None 
    for mn, mv in md.items():
        weights = np.ones_like(mv)/len(mv)
        ax = sns.distplot(mv, hist_kws={"histtype": "stepfilled"}, kde=False, rug=True, label=mn, ax=ax, norm_hist=False)

    sns.despine()
    ax = plt.axes()
    ax.set_xlabel('proportionality correlation')
    ax.set_ylabel('frequency')
    plt.legend()
    plt.savefig('{}/{}.pdf'.format(outdir, plotname))

def errorFracPlots(D, measure, methods, simIDs, outdir, plotname):#"TPRelativeError{}".format(measure))
    sns.set_style('white')
    sns.set_context('poster')
    plt.clf()
    plt.cla()

    md = {}
    for m in methods:
        md[m] = []
        for i in simIDs:
            tn = "{}_truth{}".format(measure, str(i)) 
            pn ="{}_{}{}".format(measure, m, str(i))
            tpind =  D[D[tn] >= 1]
            y = tpind[pn]
            x = tpind[tn]
            ef = 10.0
            are = 100.0 * (y - x).abs() / x
            tpef = len(are[are > ef]) / float(len(tpind))
            md[m].append(tpef)

    ax = None 
    for mn, mv in md.items():
        weights = np.ones_like(mv)/len(mv)
        ax = sns.distplot(mv, hist_kws={"histtype": "stepfilled"},kde=False, rug=True, label=mn, ax=ax, norm_hist=False)
    sns.despine()
    ax = plt.axes()
    ax.set_xlabel('TP Error Fraction')
    ax.set_ylabel('frequency')
    plt.legend()
    plt.savefig('{}/{}.pdf'.format(outdir, plotname))

   
def rsemRelDiffPlots(D, measure, methods, simIDs, outdir, plotname):
    sns.set_style('white')
    sns.set_context('poster')
    plt.clf()
    plt.cla()

    md = {}
    for m in methods:
        md[m] = []
        for i in simIDs:
            rd, _ = AnalysisUtils.relDiff("{}_truth{}".format(measure, str(i)), 
                                          "{}_{}{}".format(measure, m, str(i)), D)
            md[m].append(rd['relDiff'].abs().mean())

    ax = None 
    for mn, mv in md.items():
        weights = np.ones_like(mv)/len(mv)
        ax = sns.distplot(mv, hist_kws={"histtype": "stepfilled"},kde=False, rug=True, label=mn, ax=ax, norm_hist=False)
    sns.despine()
    ax = plt.axes()
    ax.set_xlabel('mean absolute relative difference')
    ax.set_ylabel('frequency')
    plt.legend()
    plt.savefig('{}/{}.pdf'.format(outdir, plotname))

def rsemComparisonPlots(topdir, measure, simIDs, outdir):
    import os
    kallistoQuant = "abundance.txt"
    salmonQuant = "quant.sf"
    sailfishQuant = "quant.sf/quant.sf"
    sailfishQuasiQuant = "quant.sf/quant.sf"

    print("parsing results")

    # Gather true results
    GroundTruths = []
    for i in simIDs:
        DF = pd.read_csv('{}/{}.sim.isoforms.results'.format(topdir, i), sep='\t')
        DF.drop(set(DF.columns) - set(['transcript_id', 'length', 'count', 
                'TPM', 'effective_length']), axis=1, inplace=True)
        DF.rename(columns={'transcript_id' : 'Name', 
                'length' : 'Length_truth{}'.format(i), 
                'count' : 'NumReads_truth{}'.format(i), 
                'TPM' : 'TPM_truth{}'.format(i), 
                'effective_length' : 'EffLen_truth{}'.format(i)}, inplace=True)
        DF.set_index('Name', inplace=True)
        DF.convert_objects(convert_numeric=True)
        GroundTruths.append(DF)

    SalmonRes = []
    # Gather salmon results
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'salmon', salmonQuant]) 
        DF = ParsingUtils.readSalmon(fn, '_{}{}'.format('Salmon', i))
        SalmonRes.append(DF)

    SailfishRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'sailfish', sailfishQuant]) 
        DF = ParsingUtils.readSailfish(fn, '_{}{}'.format('Sailfish', i))
        SailfishRes.append(DF)

    SalmonVBRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'salmonVB', salmonQuant])
        DF = ParsingUtils.readSalmon(fn, '_{}{}'.format('Salmon (VB)', i))
        SalmonVBRes.append(DF)

    SalmonAlnVBRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'salmon_alnVB', salmonQuant])
        DF = ParsingUtils.readSalmon(fn, '_{}{}'.format('SalmonAln (VB)', i))
        SalmonAlnVBRes.append(DF)

    SalmonAlnRes = []
    # Gather salmon results
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'salmon_aln', salmonQuant]) 
        DF = ParsingUtils.readSalmon(fn, '_{}{}'.format('SalmonAln', i))
        SalmonAlnRes.append(DF)
   
    # Gather kallisto results
    KallistoRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'kallisto', kallistoQuant]) 
        DF = ParsingUtils.readKallisto(fn, '_{}{}'.format('Kallisto', i))
        KallistoRes.append(DF)

    # Gather sailfish quasi results
    SailfishQuasiRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'sailfish_quasi', sailfishQuasiQuant]) 
        DF = ParsingUtils.readSalmon(fn, '_{}{}'.format('Sailfish (Quasi)', i))
        SailfishQuasiRes.append(DF)


    # Gather eXpress results
    expressQuant = "results.xprs"
    ExpressRes = []
    for i in simIDs:
        fn = os.path.sep.join([topdir, str(i), 'express', expressQuant])
        DF = ParsingUtils.readExpress(fn, '_{}{}'.format('eXpress', i))
        ExpressRes.append(DF)

    # Gather StringTie results
    #stringtieQuant = "t_data.ctab"
    #StringtieRes = []
    #for i in simIDs:
    #    fn = os.path.sep.join([topdir, str(i), 'stringtie', stringtieQuant])
    #    DF = ParsingUtils.readExpress(fn, '_{}{}'.format('stringtie', i))
    #    StringtieRes.append(DF)

    K = KallistoRes[0].join(KallistoRes[1:])
    S = SalmonRes[0].join(SalmonRes[1:])
    SVB = SalmonVBRes[0].join(SalmonVBRes[1:])
    SA = SalmonAlnRes[0].join(SalmonAlnRes[1:])
    SAVB = SalmonAlnVBRes[0].join(SalmonAlnVBRes[1:])
    E = ExpressRes[0].join(ExpressRes[1:])
    G = GroundTruths[0].join(GroundTruths[1:])
    SF = SailfishRes[0].join(SailfishRes[1:])
    SFQ = SailfishQuasiRes[0].join(SailfishQuasiRes[1:])
    #ST = StringtieRes[0].join(StringtieRes[1:])

    # Gather *all* results into a single dataframe
    M = G.join(K).join(S).join(SA).join(SVB).join(SAVB).join(E).join(SF).join(SFQ)

    methods = ["Salmon", "Salmon (VB)", "SalmonAln", "SalmonAln (VB)", "eXpress", "Kallisto", "Sailfish", "Sailfish (Quasi)"]

    # Filter eXpress results
    for i in simIDs:
        minVal = np.inf
        for mn in set(methods) - set(["Truth", "eXpress"]):
            print("Method name = {}".format(mn))
            newMin = M.loc[M["{}_{}{}".format(measure, mn, i)]>0, "{}_{}{}".format(measure, mn, i)].min()
            minVal = min(minVal, newMin) 
        print("filtering eXpress results < {} {}".format(minVal, measure))
        AnalysisUtils.filterValues("{}_{}{}".format(measure, "eXpress", i), M, minVal)

    #print("generating relDiff plot")
    #rsemRelDiffPlots(M, measure, methods, simIDs, outdir, "MeanAbsRelDiffs{}".format(measure))
    #print("generating proportionality plot")
    #proportionalityCorrelationPlots(M, measure, methods, simIDs, outdir, "ProportionalityCorrelation{}".format(measure))
    print("generating spearman correlation plot")
    methods_spearman = ["Sailfish", "eXpress", "Salmon",  "SalmonAln"]
    spearmanCorrelationPlots(M, measure, methods_spearman, simIDs, outdir, "SpearmanCorrelation{}".format(measure), setxlim=(0.87, 0.93))
    methods_spearman_kallisto = ["Kallisto", "Salmon"]
    spearmanCorrelationPlots(M, measure, methods_spearman_kallisto, simIDs, outdir, "SpearmanCorrelationKallisto{}".format(measure), setxlim=(0.87, 0.93))
    methods_spearman_salmon = ["Salmon (VB)","Salmon", "SalmonAln (VB)", "SalmonAln" ]
    spearmanCorrelationPlots(M, measure, methods_spearman_salmon, simIDs, outdir, "SpearmanCorrelationSalmon{}".format(measure), setxlim=(0.91, 0.93), lines=True)
    methods_spearman_sailfishquasi = ["Sailfish", "Sailfish (Quasi)", "eXpress", "Salmon",  "SalmonAln"]
    spearmanCorrelationPlots(M, measure, methods_spearman_sailfishquasi, simIDs, outdir, "SpearmanCorrelationSailfishQuasi{}".format(measure), setxlim=(0.87, 0.93), lines=True)
   #print("generating tp relative error plots")
    #errorFracPlots(M, measure, methods, simIDs, outdir, "TPErrorFrac{}".format(measure))
 
    #for i in simIDs:
    #    for m in methods + ['truth']:
    #        AnalysisUtils.filterValues("NumReads_{}{}".format(m, i), M, 1.0)
    print("generating filtered relDiff plot")
    #rsemRelDiffPlots(M, measure, methods, simIDs, outdir, "MeanAbsRelDiffs{}Filtered".format(measure))
    print("generating filtered proportionality plot")
    #proportionalityCorrelationPlots(M, measure, methods, simIDs, outdir, "ProportionalityCorrelation{}Filtered".format(measure))
    print("generating spearman correlation plot")
    #spearmanCorrelationPlots(M, measure, methods, simIDs, outdir, "SpearmanCorrelation{}Filtered".format(measure))
    #print("generating tp relative error plots")
    #tpRelativeErrorPlots(M, measure, methods, simIDs, outdir, "TPRelError{}Filtered".format(measure))
 
def plotStratifiedDiffs(M, methodDict, annotPath, outpath, measure):
    import pickle
    import seaborn as sns
    tgmap = {}
    with open('{}/tgmap.txt'.format(annotPath),'r') as ifile:
        for l in ifile:
            toks = l.split()
            tgmap[toks[0]] = toks[1]
 
    genes = pd.DataFrame([(k,v) for k,v in tgmap.items()], columns=['Name', 'Gene'])
    genes.set_index('Name', inplace=True)
    M = M.join(genes)
    
    relDiffDict = {}
    for mn, mf in methodDict.items():
        rdiffs, _ = AnalysisUtils.relDiff('{}_Truth'.format(measure), '{}_{}'.format(measure, mn), M)
        M['RelDiff_{}'.format(mn)] = rdiffs
    
    mByGenes = M.groupby('Gene')
    groups = dict(list(mByGenes))
    ntmap = {}
    methods = ["Sailfish", "Salmon", "eXpress"]
    def retainedMethod(mn):
        return mn in methods 

    sns.set_palette(sns.color_palette([methodcolor[m] for m in methods]))

    numExpressedGenes = 0
    for i,(k,group) in enumerate(groups.items()):
        if i % 10000 == 0:
            print("processing group {} of {}".format(i, len(groups)))
        ## remove 0s?
        if group['NumReads_Truth'].sum() == 0:
            continue
        numExpressedGenes += 1
        numTran = len(group)
        if numTran not in ntmap:
            ntmap[numTran] = {}
            for mn, mf in methodDict.items():
                if retainedMethod(mn):
                #if not mn.upper().endswith('TRUTH'):
                    ntmap[numTran][mn] = []
        
        for mn, mf in methodDict.items():
            if retainedMethod(mn):
            #if not mn.upper().endswith('TRUTH'):
                ntmap[numTran][mn].append(group['RelDiff_{}'.format(mn)].abs().mean())

    print("There were {} expressed genes (genes with at least 1 expressed isoform)".format(numExpressedGenes))
    nt = [set([])]
    totInClass = 0
    classSize = 1000
    for k in sorted(ntmap.keys()):
        if k > 1:
            vd = ntmap[k]
            #num = len(vd['Kallisto'])
            num = len(vd['Sailfish'])
            if len(nt[-1]) == 0 or totInClass <= classSize:
                nt[-1].add(k)
            totInClass += num
            if totInClass > classSize:
                totInClass = 0
                nt.append(set([]))
    print(nt)
    ntClasses = {}
    for c in nt:
        # make a name for the class
        minVal = min(list(c))
        maxVal = max(list(c))
        if minVal == maxVal:
            ntClasses[minVal] = str(minVal)
        else:
            name = str(minVal)+"-"+str(maxVal)
            for x in c:
                ntClasses[x] = name 
    
    from pprint import pprint
    pprint(ntClasses)

    data = []
    for k, v in ntmap.items():
        if k > 1:
            for mn, vals in v.items():
                for x in vals:
                    cname = ntClasses[k]
                    data.append((cname, mn, x))
        
    ntxp = '# txp. per gene'
    D = pd.DataFrame(data, columns=[ntxp, 'Method', 'MARDs'])
    
    sns.set_style('white')
    #sns.set_palette('pastel', desat=0.7)
    plt.clf()
    plt.cla()
    
    plt.figure(figsize=(12,8))
    cnames = sorted(list(set([v for k,v in ntClasses.items()])))
    print("Classes = ")
    pprint(cnames)
    x_order = sorted(cnames, key=lambda x : int(x.split('-')[0]))
    import matplotlib
    matplotlib.rc("lines", markersize=0, markeredgewidth=0) 
    g = sns.factorplot(ntxp, "MARDs", "Method", D, kind='box', markers='', x_order=x_order)
    #g = sns.facetgrid(ntxp, "MARDs", "Method", D)
    
    g.set_xticklabels(rotation=30)
    #g.despine(offset=10, trim=True)
    #plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('{}/GeneStratRD.pdf'.format(outpath))

#    plt.clf()
#    plt.cla()
#    
#    plt.figure(figsize=(12,8))
#    g = sns.factorplot(ntxp, "MARDs", "Method", D, kind='box', x_order=list(map(str,range(12,24))))
#    g.despine(offset=10, trim=True)
#    plt.savefig('{}/GeneStratRD{}_{}.pdf'.format(outpath, "12", "23"))
#    
#    plt.clf()
#    plt.cla()
#    
#    plt.figure(figsize=(12,8))
#    g = sns.factorplot(ntxp, "MARDs", "Method", D, kind='box', x_order=list(map(str,range(24,35)))+[u" \u226535"])
#    g.despine(offset=10, trim=True)
#    plt.savefig('{}/GeneStratRD{}_{}.pdf'.format(outpath, "24", "36"))

def makeTable(methodDict, outpath, outfile, measure, annotPath):
    import pandas as pd
    import seaborn as sns
    import ParsingUtils
    import AnalysisUtils
    dframes = []
    for k, v in methodDict.items():
        if k.upper().startswith('SALMON'):
            d = ParsingUtils.readSalmon(v, '_{}'.format(k))
        elif k.upper().startswith('KALLISTO'):
            d = ParsingUtils.readKallisto(v, '_{}'.format(k))
        elif k.upper().startswith('EXPRESS'):
            d = ParsingUtils.readExpress(v, '_{}'.format(k))
        elif k.upper() == 'SAILFISH':
            d = ParsingUtils.readSailfish(v, '_{}'.format(k))
        elif k.upper() == 'SAILFISH (QUASI)':
            d = ParsingUtils.readSalmon(v, '_{}'.format(k))
        elif k.upper().startswith('TRUTH'):
            suffix = '_{}'.format(k)
            d = ParsingUtils.readProFile(v, suffix) 
            d["TPM{}".format(suffix)] = 1000000.0 * (d["ExpFrac{}".format(suffix)] / d["ExpFrac{}".format(suffix)].sum())
            # Flux sim thinks paired-end = 2 reads . . . sinh
            d["NumReads{}".format(suffix)] = d["SeqNum{}".format(suffix)] * 0.5

        # Add this dataframe to the list
        dframes.append(d)

    M = dframes[0].join(dframes[1:])
    
    # Filter eXpress results
    minVal = np.inf
    for mn in set(methodDict.keys()) - set(["Truth", "eXpress"]):
        newMin = M.loc[M["{}_{}".format(measure, mn)]>0, "{}_{}".format(measure,mn)].min()
        minVal = min(minVal, newMin) 
    print("filtering eXpress results < {} {}".format(minVal, measure))
    AnalysisUtils.filterValues("{}_{}".format(measure, "eXpress"), M, minVal)

    org = outfile.split('/')[-1].split('_')[0] 
    print("org = {}".format(org))
    if org == 'human':
        plotStratifiedDiffs(M, methodDict, annotPath, outpath, measure)

    mrdName = 'abs. mean rel. diff.'
    corrName = 'Spearman corr.'
    propName = 'Proportionality corr.'
    tpefName = 'TP error fraction'
    tpMedErrorName = 'TP median per. error'
    res = pd.DataFrame(data={ m : {tpMedErrorName : np.nan, tpefName : np.nan, mrdName : np.nan, corrName : np.nan, propName : np.nan} for m in (methodDict.keys() - set('Truth'))})

    import scipy as sp
    import scipy.stats

    for k in methodDict:
        if k.upper() != "TRUTH":
            c = sp.stats.spearmanr(M["{}_Truth".format(measure)], M["{}_{}".format(measure, k)])[0]
            res[k][corrName] = c
            mrd, _ = AnalysisUtils.relDiff("{}_Truth".format(measure), "{}_{}".format(measure, k), M) 
            res[k][mrdName] = mrd["relDiff"].abs().mean()

            pc = AnalysisUtils.proportionalityCorrelation("{}_Truth".format(measure), "{}_{}".format(measure, k), M) 
            res[k][propName] = pc 

            tpind =  M[M["{}_Truth".format(measure)] >= 1]
            y = tpind["{}_{}".format(measure, k)] 
            x = tpind["{}_Truth".format(measure)]
            ef = 10.0
            re = (y - x) / x
            are = 100.0 * (y - x).abs() / x
            tpef = len(are[are > ef]) / float(len(are))
            res[k][tpefName] = tpef
            res[k][tpMedErrorName] = re.median()

    res.drop('Truth', axis=1, inplace=True)
    print(res)
    res.to_csv(outfile+".csv")

    with open(outfile, 'w') as ofile:
        ofile.write(res.to_latex(float_format=lambda x: "{0:.2f}".format(x)))
    print("wrote {}".format(outpath))

def true_neg(tc, pc, d, cutoff=1e-2):
    return len(d.loc[d[tc] < cutoff][d[pc] < cutoff])

def true_pos(tc, pc, d, cutoff=1e-2):
    return len(d.loc[d[tc] >= cutoff][d[pc] >= cutoff])

def false_neg(tc, pc, d, cutoff=1e-2):
    return len(d.loc[d[tc] >= cutoff][d[pc] < cutoff])

def false_pos(tc, pc, d, cutoff=1e-2):
    return len(d.loc[d[tc] < cutoff][d[pc] >= cutoff])




def makePlots(truthFile, predFile, method, outBase, measure):

    if measure == "tpm":
        trueColName = "TPM_truth"
        predColName = "TPM_{}".format(method)
    elif measure == "num_reads":
        trueColName = "NumReads_truth"
        predColName = "NumReads_{}".format(method)

    # Load the data; first the predicitons
    p =  None
    if method == "salmon":
        p = ParsingUtils.readSalmon(predFile, '_{}'.format(method))
    elif method == "kallisto":
        p = ParsingUtils.readKallisto(predFile, '_{}'.format(method))
    elif method == "express":
        p = ParsingUtils.readExpress(predFile, '_{}'.format(method))
    
    # Now the ground truth
    g = ParsingUtils.readProFile(truthFile, '_truth')

    # Convert to TPM
    g["TPM_truth"] = 1000000.0 * (g["ExpFrac_truth"] / g["ExpFrac_truth"].sum())
    # Flux sim thinks paired-end = 2 reads . . . sigh
    g["NumReads_truth"] = g["SeqNum_truth"] * 0.5

    # Filter out low TPM
    AnalysisUtils.filterValues("TPM_truth", g, 0.01)
    AnalysisUtils.filterValues("TPM_{}".format(method), p, 0.01)
    AnalysisUtils.filterValues("NumReads_truth", g, 1.0)
    AnalysisUtils.filterValues("NumReads_{}".format(method), p, 1.0)

    # merge dataframes
    m = g.join(p)

    setPlotProperties()

    makeCorrPlot(trueColName, predColName, m, outBase, method, measure)
    makeRelDiffPlot(trueColName, predColName, m, outBase, method, measure)
    makeRelErrorPlot(trueColName, predColName, m, outBase, method, measure)
    
def get_time(fname):
    import json
    with open(fname) as f:
        return json.load(f)['wall_clock_times']['s'][0]

def prettify_org(org):
    human_pretty = 'Human'
    zeamays_pretty = 'Zea mays'
    if org == 'human':
        org = human_pretty
    elif org == 'zeamays':
        org = zeamays_pretty
    return org

def prettify_method(method):
    if len(method) == 1:
        if method[0] == "express": 
            method = "eXpress"
        else:
            method = method[0].capitalize()
    else:
        print("method = {}".format(method))
        suffix = ''
        if method[0] == "useVBOpt": suffix = " (VB)"
        elif method[-1] == "quasi": suffix = " (Quasi)"

        if suffix == " (Quasi)":
            method = ''.join([x.capitalize() for x in method[:-1]]) + suffix
        else:
            method = ''.join([x.capitalize() for x in method[1:]]) + suffix

    return method

# Different naming convention for RSEM (maybe change later)
def prettify_method_rsem(method):
    print("method[0] = {}".format(method[0]))
    if method[0] == "star":
        return "STAR Alignment"
    elif method[0] == "express": 
        return "eXpress"
    else:
        method[0] = method[0].capitalize()
        if method[0].find("vb") != -1:
            method[0] = method[0].rstrip("vb") + " (VB)"

    if len(method) > 1:
        if method[1].find("aln") != -1:
            method[0] += "Aln"
        if method[1].find("VB") != -1:
            method[0] += " (VB)"

    print("method = {}".format(method[0]))
    return method[0]

def makeBenchmarkPlot(timing_files, timing_files_alignments, outfile):
    import glob, json
    font = {
            'weight' : 'bold',
            'size'   : 50}

    import matplotlib
    matplotlib.rc('font', **font)

    sns.set_style("whitegrid")
    plt.cla()
    plt.clf()

    methodDict = {}
    timingData = []


    for fname in timing_files:
        numsecs = get_time(fname)
        pieces = fname.split('.')[0].split("/")[-1].split("_")
        org = pieces[0]
        method = pieces[1:]
        org = prettify_org(org)
        method = prettify_method(method)

        timingData.append( (org, method, numsecs) )

    for fname in timing_files_alignments:
        numsecs = get_time(fname)
        if fname.find("human") != -1:
            org = prettify_org("human")
        elif fname.find("zeamays") != -1:
            org = prettify_org("zeamays")

        method = "STAR Alignment"
        timingData.append( (org, method, numsecs) )

    #x_order = ['Kallisto', 'Salmon', 'Salmon (VB)', 'SalmonAln', 'SalmonAln (VB)', 'STAR Alignment', 'eXpress']
    #x_order = ['Sailfish', 'Salmon', 'STAR Alignment', 'eXpress']
    x_order = ['Sailfish', 'Salmon', 'STAR Alignment', 'SalmonAln']
    benchmarkDF = pd.DataFrame(timingData, columns=['Organism', 'Method', 'Time (s)']) 
    print( benchmarkDF ) 
    plt.figure(figsize=(10,6))
    sns.barplot("Method", "Time (s)", "Organism", data=benchmarkDF, x_order=x_order)
    sns.despine()
    plt.savefig(str(outfile))    

def makeBenchmarkPlotSupp(timing_files, outfile):
    import glob, json
    font = {
            'weight' : 'bold',
            'size'   : 50}

    import matplotlib
    matplotlib.rc('font', **font)

    sns.set_style("whitegrid")
    plt.cla()
    plt.clf()

    methodDict = {}
    timingData = []


    for fname in timing_files:
        numsecs = get_time(fname)
        pieces = fname.split('.')[0].split("/")[-1].split("_")
        org = pieces[0]
        method = pieces[1:]
        org = prettify_org(org)
        method = prettify_method(method)
        print(method)
        timingData.append( (org, method, numsecs) )

    #x_order = ['Kallisto', 'Salmon', 'Salmon (VB)', 'SalmonAln', 'SalmonAln (VB)', 'STAR Alignment', 'eXpress']
    x_order = ['Sailfish', 'Sailfish (Quasi)', 'Salmon']#'Salmon', 'STAR Alignment', 'eXpress']
    benchmarkDF = pd.DataFrame(timingData, columns=['Organism', 'Method', 'Time (s)']) 
    print( benchmarkDF ) 
    plt.figure(figsize=(10,6))
    sns.barplot("Method", "Time (s)", "Organism", data=benchmarkDF, x_order=x_order)
    sns.despine()
    plt.savefig(str(outfile))    

def makeBenchmarkPlotSupp2(timing_files, outfile):
    import glob, json
    font = {
            'weight' : 'bold',
            'size'   : 50}

    import matplotlib
    matplotlib.rc('font', **font)

    sns.set_style("whitegrid")
    plt.cla()
    plt.clf()

    methodDict = {}
    timingData = []


    for fname in timing_files:
        numsecs = get_time(fname)
        pieces = fname.split('.')[0].split("/")[-1].split("_")
        org = pieces[0]
        method = pieces[1:]
        org = prettify_org(org)
        method = prettify_method(method)
        print(method)
        timingData.append( (org, method, numsecs) )

    #x_order = ['Kallisto', 'Salmon', 'Salmon (VB)', 'SalmonAln', 'SalmonAln (VB)', 'STAR Alignment', 'eXpress']
    x_order = ['Kallisto', 'Salmon']#'Salmon', 'STAR Alignment', 'eXpress']
    benchmarkDF = pd.DataFrame(timingData, columns=['Organism', 'Method', 'Time (s)']) 
    print( benchmarkDF ) 
    plt.figure(figsize=(10,6))
    sns.barplot("Method", "Time (s)", "Organism", data=benchmarkDF, x_order=x_order)
    sns.despine()
    plt.savefig(str(outfile))    


def makeRSEMBenchmarkPlot(replicates, outfile):
    import glob
    sns.set_style("whitegrid")
    plt.cla()
    plt.clf()
    timingData = {}
    dataList = []
    for r in replicates:
        method_fnames = glob.glob("{}/*.json".format(r))
        for fname in method_fnames:
            numsecs = get_time(fname)
            method = fname.split("/")[-1].split(".")[0].split("_")
            method = prettify_method_rsem(method)
            dataList.append( (method, numsecs) )
            if method not in timingData:
                timingData[method] = []
            timingData[method].append(numsecs)
    benchmarkDF = pd.DataFrame(dataList, columns=['Method', 'Time (s)'])
    x_order = ['Kallisto', 'Salmon', 'Salmon (VB)', 'SalmonAln', 'SalmonAln (VB)', 'STAR Alignment', 'eXpress']
    plt.figure(figsize=(10,6))
    sns.barplot("Method", "Time (s)", data=benchmarkDF, x_order=x_order)
    ax = plt.axes()
    sns.despine()
    plt.savefig(str(outfile))


def qpcrplots(outpath):
    font = {
            'weight' : 'bold',
            'size'   : 50}
    import matplotlib
    matplotlib.rc('font', **font)
    sns.set_style("white")
    tgmap = pd.read_csv('/shared/SalmonManuscript/data/HumanGenomeAnnotations/tgmap.all.txt', sep='\t', names=['gene.id', 'Name', 'gene.name'])
    tgmap.set_index('Name', inplace=True)
    print( "length tgmap={}".format(len(tgmap)) )
    salmondf = ParsingUtils.readSalmon('data/quant_seqc/rep1_salmon/quant.sf', suffix="_salmon")
    print( "length salmondf={}".format( len(salmondf)) )
    merged = salmondf.join(tgmap, how='inner')
    print( "length tgmap&salmondf={}".format( len(merged) ) )
    expressDF = ParsingUtils.readExpress('data/quant_seqc/rep1_express/results.xprs', suffix="_express")
    merged = merged.join(expressDF, how='inner')
    print( "Length tgmap&salmonDF&expressDF={}".format( len(merged) ) )
    sailfishDF = ParsingUtils.readSailfish('data/quant_seqc/rep1_sailfish/quant_bias_corrected.sf', suffix="_sailfish")
    merged = merged.join(sailfishDF, how='inner')
    print( "Length tgmap&salmonDF&expressDF&sailifishDF={}".format( len(merged) ) )
    kallistoDF = ParsingUtils.readKallisto('data/quant_seqc/rep1_kallisto/abundance.txt', suffix="_kallisto")
    merged = merged.join(kallistoDF, how='inner')
    print( "Length tgmap&salmonDF&expressDF&sailifishDF&kallistoDF={}".format( len(merged) ) )
    merged = merged.groupby('gene.name').sum()
    print( "Length after summing by gene={}".format(len(merged)) ) 
    qpcrDF = pd.read_csv('data/seqc/qpcr/PrimePCR.txt', sep='\t', names=['gene.name', 'A','B','C','D'], skiprows=1)
    qpcrDF.set_index("gene.name", inplace=True)
    print( "Length qpcrDF={}".format(len(qpcrDF)) )
    merged = merged.join(qpcrDF, how='inner')
    print( "Length everything&qpcrDF={}".format(len(merged)) )

    methods = ["Sailfish", "Kallisto", "Salmon", "eXpress"]
    #sns.set_palette(sns.color_palette([methodcolor[m] for m in methods]))

    import AnalysisUtils
    from AnalysisUtils import proportionalityCorrelation 
    prop_salmon = proportionalityCorrelation('TPM_salmon', 'A', merged)
    prop_express = proportionalityCorrelation('TPM_express', 'A', merged)
    prop_sailfish = proportionalityCorrelation('TPM_sailfish', 'A', merged)
    prop_kallisto = proportionalityCorrelation('TPM_kallisto', 'A', merged)
    print( "Proportionality Salmon: {}".format(prop_salmon) )
    print( "Proportionality eXpress: {}".format(prop_express) )
    print( "Proportionality Sailfish: {}".format(prop_sailfish) )
    print( "Proportionality Kallisto: {}".format(prop_kallisto) )

    spearman_salmon =  merged['TPM_salmon'].corr( merged['A'], "spearman")
    spearman_express =  merged['TPM_express'].corr( merged['A'], "spearman")
    spearman_sailfish =  merged['TPM_sailfish'].corr( merged['A'], "spearman")
    spearman_kallisto =  merged['TPM_kallisto'].corr( merged['A'], "spearman")
    print( "Spearman Salmon: {}".format(spearman_salmon) )
    print( "Spearman eXpress: {}".format(spearman_express) )
    print( "Spearman Sailfish: {}".format(spearman_sailfish) )
    print( "Spearman Kallisto: {}".format(spearman_kallisto) )
    #plt.xlim(0,17000)
    def plot_table(data, methods, suffix=""):
        plt.clf()
        plt.figure(figsize=(8,12))
        corrtable = pd.DataFrame(data=data, columns = ['Spearman corr.'], index=methods)
        corrtable = corrtable.reset_index()
        corrtable.columns = ['Method', 'Spearman corr.']
        print (corrtable)
        sns.barplot("Method", "Spearman corr.", data=corrtable)
        plt.savefig("{}/qpcrbarplot{}.pdf".format(outpath,suffix))
    methods = ["Sailfish", "Salmon", "eXpress"]
    data = [spearman_sailfish, spearman_salmon, spearman_express]
    plot_table(data,methods)
    methods = ["Kallisto", "Salmon"]
    data = [spearman_kallisto, spearman_salmon]
    plot_table(data,methods,"_supp")
    plt.cla()
    plt.clf()
    plt.plot(merged['TPM_salmon'].rank()[7000:14000], merged['A'].rank()[7000:14000], 'o', alpha=1, ms=2)
    plt.savefig("{}/qpcr_salmon.pdf".format(outpath))
    plt.clf()
    plt.plot(merged['TPM_express'].rank()[7000:14000], merged['A'].rank()[7000:14000], 'o',alpha=1, ms=2)
    #plt.xlim(0,17000)
    plt.savefig("{}/qpcr_express.pdf".format(outpath))
    plt.clf()
    plt.plot(merged['TPM_sailfish'].rank()[7000:14000], merged['A'].rank()[7000:14000], 'o', alpha=1, ms=2)
    #plt.xlim(0,17000)
    plt.savefig("{}/qpcr_sailfish.pdf".format(outpath))
    plt.clf()
    plt.plot(merged['TPM_kallisto'].rank()[7000:14000], merged['A'].rank()[7000:14000], 'o', alpha=1, ms=2)
    #plt.xlim(0,17000)
    plt.savefig("{}/qpcr_kallisto.pdf".format(outpath))

