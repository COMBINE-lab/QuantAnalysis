# Sailfish Analysis

A collection of analyses of the results of Sailfish (and other quantifiers) on a growing collection of synthetic and real data.

## An analysis in simulated data from D. *mel* from [Soneson et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3) paper.

For a further description of how these are computed, and a deeper analysis, take a look at the related [Jupyter notebook](https://github.com/COMBINE-lab/QuantAnalysis/blob/master/analysis_scripts/AnalyzeSoneson.ipynb).

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
using Sailfish `0.9.1` and kallisto `0.42.5`
