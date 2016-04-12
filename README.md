# Sailfish Analysis

A collection of analyses of the results of Sailfish (and other quantifiers) on a growing collection of synthetic and real data.

## An analysis in simulated data from the [Soneson et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3) paper.


### *D.* *mel* 

For a further description of how these are computed, and a deeper analysis, take a look at the related [Jupyter notebook](https://github.com/COMBINE-lab/QuantAnalysis/blob/master/analysis_scripts/AnalyzeSonesonDmel.ipynb).

```
mean of medians of signed relative differences is :
 kallisto: -0.04
 sailfish: 0.00

median of means of signed relative differences is :
 kallisto: -0.04
 sailfish: 0.00

mean of medians of abs relative differences is :
 kallisto: 0.09
 sailfish: 0.06

median of means of abs relative differences is :
 kallisto: 0.29
 sailfish: 0.28
```

### *H.* *sapiens*

For a further description of how these are computed, and a deeper analysis, take a look at the related [Jupyter notebook](https://github.com/COMBINE-lab/QuantAnalysis/blob/master/analysis_scripts/AnalyzeSonesonHs.ipynb).

```
kallisto_pearson: 1.0,0.99,1.0,1.0,1.0,1.0
kallisto_spearman: 0.92,0.92,0.92,0.92,0.92,0.92

sailfish_pearson: 1.0,1.0,1.0,1.0,1.0,1.0
sailfish_spearman: 0.92,0.92,0.92,0.92,0.92,0.92

mean of medians of signed relative differences is :
 kallisto: 0.00
 sailfish: 0.00

median of means of signed relative differences is :
 kallisto: -0.01
 sailfish: -0.01

mean of medians of abs relative differences is :
 kallisto: 0.00
 sailfish: 0.00

median of means of abs relative differences is :
 kallisto: 0.23
 sailfish: 0.23
```

using Sailfish `0.10.0` and kallisto `0.42.5`


## A re-analysis of the RSEM simulated data from the [kallisto](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3519.html) paper

For a further description of how these are computed, and a deeper analysis, take a look at the related [Jupyter notebook](https://github.com/COMBINE-lab/QuantAnalysis/blob/master/analysis_scripts/AnalyzeRSEM.ipynb).

```
mean of medians of signed relative differences is :
 kallisto: 0.00
 sailfish: 0.00

median of means of signed relative differences is :
 kallisto: 0.02
 sailfish: 0.02

mean of medians of abs relative differences is :
 kallisto: 0.01
 sailfish: 0.00

median of means of abs relative differences is :
 kallisto: 0.28
 sailfish: 0.28
```
using Sailfish `0.10.0` and kallisto `0.42.5`
