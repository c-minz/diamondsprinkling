MATLAB Source Code for Simulations
==================================

This repository is part of my PhD project at the Department of Mathematics, University of York, from 2017 to 2021. The project was funded by the Engineering and Physical Sciences Research Council (EPSRC). The simulations in this project were undertaken on the Viking Cluster, which is a high performance compute facility provided by the University of York.

The results of this project are available at https://arxiv.org/abs/2011.02965

Project Information
===================

This part of my PhD project investigates the local structure in causal sets that are sprinkled in Alexandrov subsets of Minkowski spacetime. By local structure, we mean a preferred past structure (an additional structure to discretise the Klein-Gordon field equations for causal set theory) [1], and diamonds (intervals from events to elements in their rank 2 past that can be chosen as a preferred past). The aim of the simulation part of the project was to examine methods to obtain a preferred past structure from a given causal set without manually imposing it. Find more details about the research in our publication, https://arxiv.org/abs/2011.02965 .

[1] E. Dable-Heath, C. J. Fewster, K. Rejzner, and N. Woods, Algebraic classical and quantum field theory on causal sets, Phys. Rev. D 101, 065013 (2020). https://doi.org/10.1103/PhysRevD.101.065013

Structure of the Source Code
============================

The research results are generated with the following steps:
1. Raw data files are generated with the MATLAB functions and scripts provided in _this repository_. The function "runjobs\_generate\_statistics.m" is called via a job script like "N6000D*.job".
2. Raw data files are merged to a single data file (per parameter setting) with the MATLAB functions and scripts provided in _this repository_. The function "runjobs\_combine\_statistics.m" is called via a job script like "cbN6000D*.job".
3. Raw data files are converted to data tables (two column csv files) with the functions given in the repository https://github.com/c-minz/diamondresults .
4. Data files are plotted with histograms and diagrams with TikZ, also see repository https://github.com/c-minz/diamondresults .

Distribution Information
========================

Development environments: MATLAB R2019a, MATLAB R2020a

Software version on the cluster: MATLAB R2018a

Copyright 2021 by Christoph Minz.