# pgsbenchmark.org

The platform for privacy-preserving benchmarks for polygenic prediction.
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev) -->
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev?labpath=minimal.ipynb) -->

## Leaderboard



## Abstract of accompanying paper

Recently, several new approaches for creating polygenic scores (PGS) have been developed and this trend shows no sign of abating. However, it has thus far been challenging to determine which approaches are superior, as different studies report seemingly conflicting benchmark results. This heterogeneity in benchmark results is in part due to different outcomes being used, but also due to differences in the genetic variants being used, data preprocessing, and other quality control steps. As a solution, a publicly available benchmark for polygenic prediction is presented here, which allows researchers to both **train** and **test** polygenic prediction methods using only summary-level information, thus preserving privacy. Using simulations and real data, we show that model performance can be estimated with accuracy, using only linkage disequilibrium (LD) information and genome-wide association summary statistics for target outcomes. Finally, we make this PGS benchmark - consisting of 8 outcomes, including somatic and psychiatric disorders - publicly available for researchers to download on our PGS benchmark platform (http://www.pgsbenchmark.org). We believe this benchmark can help establish a clear and unbiased standard for future polygenic score methods to compare against.

In other words, the benchmark data on this site provides everything needed to develop polygenic score methods yourself. We hope that this resource will also become a great resource for determining which PGS approaches are state-of-the-art if many authors of PGS approaches use the resource to evaluate their methods and email us the result.

Link to the paper:
https://www.biorxiv.org/content/10.1101/2022.10.10.510645v1

<br><br>

## Condensed results

![image](https://user-images.githubusercontent.com/6292714/195577590-a8b9e900-bcd8-41ae-a7d6-edc42322cb35.png)


## The Benchmark data

The files contain all the data that is needed for computing the benchmark and develop PGS methods yourself.
It is save in HDF5 format.

External Sumstats & Validation dataset (LD + matched GWAS sumstats for the 10K induviduals):
https://drive.google.com/file/d/1eZIiuz__tiqEtTyL14C_GwPEC6raSFSf

Test dataset (LD + matched GWAS sumstats for 352K induviduals):
https://drive.google.com/file/d/12v67vJoAZqkvaStqiYzUrAezS4KTNdNP

I will add a tutorial with code soon to make using the data easier.


