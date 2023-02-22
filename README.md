# pgsbenchmark.org

<!-- The best way to look at this information is to go to pgsbenchmark.org!
Some additional information for the "Code and Software Submission Checklist" can be found
at the end of this document in a html comment like this one. -->

The platform for privacy-preserving benchmarks for polygenic prediction. 

It is potentially useful for: 
1.	**People that develop PGS methods ->** The data on pgsbenchmark.org allows them to benchmark their approach and provides everything needed for developing a PGS method.
2.	**People that are looking for a PGS method to use ->** They can see which method is currently state-of-the-art.

<br>

<!-- notes to self:
[Do_we_even_need_UCB.ipynb](https://colab.research.google.com/github/SMPyBandits/SMPyBandits/blob/master/notebooks/Do_we_even_need_UCB.ipynb) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/dev?labpath=nbs/loaders.ipynb)
[![Google Colab](https://badgen.net/badge/Launch/on%20Google%20Colab/blue?icon=terminal)](https://colab.research.google.com/github/SMPyBandits/SMPyBandits/blob/master/notebooks/Do_we_even_need_UCB.ipynb)
https://github.com/mennowitteveen/pgsbenchmark/blob/dev/nbs/loaders.ipynb
https://github.com/binder-examples/jupyter-extension
-->


## The Leaderboard

This table gives the result for the benchmark and tells which approaches are state-of-the-art. More approaches will be added later.
If you want to have your method on the leaderboard, please send me an email.

<img width="218" alt="image" src="https://user-images.githubusercontent.com/6292714/198376342-38f04f6f-c2cd-481a-8072-37b11e9c0625.png">
<br> 

## Getting started

Run the following line on a Unix system with Jupyter and Python (3.2+) installed.
```
git clone https://github.com/mennowitteveen/pgsbenchmark.git
pip install -r ./pgsbenchmark/requirements.txt
jupyter notebook ./pgsbenchmark/nbs/PPB-demonstration.ipynb
```

or do it the easy way with Google Colab, by clicking the following link:

[![Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mennowitteveen/pgsbenchmark/blob/main/nbs/PPB-demonstration.ipynb)

or if the Colab messages *"Warning: This notebook was not authored by Google."* scared you, you can alternatively click the following button for a similar cloud experience (=slower though):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mennowitteveen/pgsbenchmark/main?labpath=nbs/PPB-demonstration.ipynb)

<br>

## The Benchmark data

The files contain all the data that is needed for computing the benchmark and develop PGS methods yourself.
It is saved in HDF5 format. 

**External Sumstats & Validation dataset (LD + matched GWAS sumstats for the 10K induviduals):**<br>
https://drive.google.com/file/d/1eZIiuz__tiqEtTyL14C_GwPEC6raSFSf

**Test dataset (LD + matched GWAS sumstats for 352K induviduals):**<br>
https://drive.google.com/file/d/12v67vJoAZqkvaStqiYzUrAezS4KTNdNP

<br>

## Abstract of accompanying paper

Recently, several new approaches for creating polygenic scores (PGS) have been developed and this trend shows no sign of abating. However, it has thus far been challenging to determine which approaches are superior, as different studies report seemingly conflicting benchmark results. This heterogeneity in benchmark results is in part due to different outcomes being used, but also due to differences in the genetic variants being used, data preprocessing, and other quality control steps. As a solution, a publicly available benchmark for polygenic prediction is presented here, which allows researchers to both **train** and **test** polygenic prediction methods using only summary-level information, thus preserving privacy. Using simulations and real data, we show that model performance can be estimated with accuracy, using only linkage disequilibrium (LD) information and genome-wide association summary statistics for target outcomes. Finally, we make this PGS benchmark - consisting of 8 outcomes, including somatic and psychiatric disorders - publicly available for researchers to download on our PGS benchmark platform (http://www.pgsbenchmark.org). We believe this benchmark can help establish a clear and unbiased standard for future polygenic score methods to compare against.

In other words, the benchmark data on this site provides everything needed to develop polygenic score methods yourself. We hope that this resource will also become a great resource for determining which PGS approaches are state-of-the-art if many authors of PGS approaches use the resource to evaluate their methods and email us the result.


**Witteveen, M.J., Pedersen, E.M., Meijsen, J., Andersen, M.R., Privé, F., Speed, D. and Vilhjálmsson, B.J., 2022. Publicly Available Privacy-preserving Benchmarks for Polygenic Prediction. bioRxiv.**


Link to the paper:
https://www.biorxiv.org/content/10.1101/2022.10.10.510645v1

<br>

![image](https://user-images.githubusercontent.com/6292714/195577590-a8b9e900-bcd8-41ae-a7d6-edc42322cb35.png)


<br>
Some extra information that will not be needed by most users can be found in the raw text of this readme.
<br>


<!--
Extended Readme information:

System requirements:

- The all this code should be run using python (3.2+) on a Unix system (Linux or Mac). Required Python packages can be found in requirements.txt
- The versions that the software has actually been tested on can be found in usedversions.txt
- There are no special hardware requirements for running the code.

Installation guide:

- Installation instructions are provided above in the section "Getting started".
- Installation time on a normal machine will be less then 2 minutes.

Demo:

- To run the demo execute the notebook nbs/PPB-demonstration.ipynb using Jupyter.
- The expected output is what is already present in the notebook. This should reproduce exactly.
- The expected runtime is a couple of mi nutes on most systems.

Instructions for use:

- The Privacy preserving benchmark data can be loaded and used as in the demonstration.
- The other code in nbs/ was used to generate the results in the paper.
- Full replication of the paper results requires individual-level data to be included (UKBB+iPSYCH) which could not be included, since this genetic data is restricted.

-->



