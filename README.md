# Evolutionary history of the *Anolis distichus* species group

This repository contains the code necessary to reproduce the results of Myers et al. *Adaptive signal evolution and gene flow drive variable
speciation outcomes in bark anoles*.

## Bioinformatics

| ddRADseq filtering and population clustering visualization | [scripts/Sampling_PCA.Rmd](scripts/Sampling_PCA.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/Sampling_PCA.html) |

## Population genetics

| Population differentiation (F<sub>ST</sub>) & diversity (π) estimation and IBD analysis  | [scripts/Pairwise-Fst-IBD.Rmd](scripts/Pairwise-Fst-IBD.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/Pairwise-Fst-IBD.html) |

| Population network estimation with NeighborNet algorithm with `phangorn`  | [scripts/neighborNet.Rmd](scripts/neighborNet.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/neighborNet.html) |

## Phylogenetics

| Concatenated ddRAD loci maximum likelihood phylogeny | [scripts/concatenated-ML-tree.rmd](scripts/concatenated-ML-tree.rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/concatenated-ML-tree.html) |

| Mitochondrial DNA phylogeny | [scripts/mtdna.Rmd](scripts/mtdna.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/mtdna.html) |

| Sex chromosome--autosome discordance | [scripts/y_chromosome_analyses.Rmd](scripts/y_chromosome_analyses.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/y_chromosome_analyses.html) |

| Species tree inference with `SVDquartets` | [scripts/svdquartets.rmd](scripts/svdquartets.rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/svdquartets.html) |

| Species tree inference using polymorphism-aware phylogenetic model (PoMo) | [scripts/PoMo-species-tree-estimation.rmd](scripts/PoMo-species-tree-estimation.rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/PoMo-species-tree-estimation.html) |


## Gene flow Analyses

| Allele-sharing patterns with `Dsuite`  | [scripts/Dsuite.Rmd](scripts/Dsuite.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/Dsuite.html) |

| Estimating migration rates with the MSC-M in `BPP`  | [scripts/migration.Rmd](scripts/migration.Rmd) | [HTML](https://TannerMyers.github.io/bark-anole-evolution/migration.html) |
