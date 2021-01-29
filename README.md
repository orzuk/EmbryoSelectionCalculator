# EmbryoSelectionCalculator
Calculator for the potential benefit of embryo selection experiments using Polygenic Risk Scores


This repository contains R code used to calculate the possible gain in phenotype value achievable in embryo selection experiments 
under Gaussian models for the genetic architecture. 
The repository accompanies the paper [1]. Please cite this paper if you're using the package. 

## Installation
clone the repository into a directory of your choice, and start an R session within this directory. 

## Usage example 
G <- gain.moments(10, 0.3, 'exact', 1)

G$E.G will contain the expected gain (in units of standard deviations) in selecting the top quantitative trait out of n=10 embryos, for selection based on a polygenic score which explaines r2ps=0.3 of the trait's variance. 
G$Var.G will contain the variance of this gain. 


### Authors
For any questions, please contact Shai Carmi (shai.carmi@mail.huji.ac.il) or Or Zuk (or.zuk@mail.huji.ac.il)


### Ref
[1] Screening Human Embryos for Polygenic Traits has Limited Utility <br>
E. Karavani, O. Zuk, D. Zeevi, G. Atzmon, N. Barzilai, N.C. Stefanis, A. Hatzimanolis, N. Smyrnis, D. Avramopoulos, L. Kruglyak, M. Lam, T. Lencz and S. Carmi <br>
https://www.sciencedirect.com/science/article/pii/S0092867419312103  , Cell 179(6), p. 1424-1435 (2019) 
 
