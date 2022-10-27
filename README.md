# pgsbenchmark.org

The platform for privacy-preserving benchmarks for polygenic prediction. ]

It is potentially useful for: 
1.	**People that develop PGS methods -> t**he data on pgsbenchmark.org allows them to benchmark their approach and provides everything needed for developing a PGS method.
2.	**People that are looking for a PGS method to use ->** They can see which method is currently state-of-the-art (especially if point 1 happens).

<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev) -->
<!-- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev?labpath=minimal.ipynb) -->

## The Leaderboard

This table gives the result for the benchmark and tells which approaches are state-of-the-art. More approaches will be added later.
If you want to have your method on the leaderboard, please send me an email.

<img width="218" alt="image" src="https://user-images.githubusercontent.com/6292714/198376342-38f04f6f-c2cd-481a-8072-37b11e9c0625.png">
<br>

## The Benchmark data

The files contain all the data that is needed for computing the benchmark and develop PGS methods yourself.
It is save in HDF5 format.

**External Sumstats & Validation dataset (LD + matched GWAS sumstats for the 10K induviduals):**<br>
https://drive.google.com/file/d/1eZIiuz__tiqEtTyL14C_GwPEC6raSFSf

**Test dataset (LD + matched GWAS sumstats for 352K induviduals):**<br>
https://drive.google.com/file/d/12v67vJoAZqkvaStqiYzUrAezS4KTNdNP

I will add a tutorial with code soon to make using the data easier.

<br>

## Abstract of accompanying paper

Recently, several new approaches for creating polygenic scores (PGS) have been developed and this trend shows no sign of abating. However, it has thus far been challenging to determine which approaches are superior, as different studies report seemingly conflicting benchmark results. This heterogeneity in benchmark results is in part due to different outcomes being used, but also due to differences in the genetic variants being used, data preprocessing, and other quality control steps. As a solution, a publicly available benchmark for polygenic prediction is presented here, which allows researchers to both **train** and **test** polygenic prediction methods using only summary-level information, thus preserving privacy. Using simulations and real data, we show that model performance can be estimated with accuracy, using only linkage disequilibrium (LD) information and genome-wide association summary statistics for target outcomes. Finally, we make this PGS benchmark - consisting of 8 outcomes, including somatic and psychiatric disorders - publicly available for researchers to download on our PGS benchmark platform (http://www.pgsbenchmark.org). We believe this benchmark can help establish a clear and unbiased standard for future polygenic score methods to compare against.

In other words, the benchmark data on this site provides everything needed to develop polygenic score methods yourself. We hope that this resource will also become a great resource for determining which PGS approaches are state-of-the-art if many authors of PGS approaches use the resource to evaluate their methods and email us the result.

Link to the paper:
https://www.biorxiv.org/content/10.1101/2022.10.10.510645v1

<br>

![image](https://user-images.githubusercontent.com/6292714/195577590-a8b9e900-bcd8-41ae-a7d6-edc42322cb35.png)





