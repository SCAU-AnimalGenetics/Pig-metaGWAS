# Analysis pipelines for meta-GWAS

## 1. Introduction

As an essential part of the ongoing Farm animal Genotype-Tissue Expression (FarmGTEx) - PigGTEx (http://piggtex.farmgtex.org/) project, the meta-GWAS aims to build a a comprehensive public catalogue of genetic association map for complex traits in pigs from diverse regional distribution. This resources will be updated every 2-3 years together with PigGTEx project. The meta-GWAS will serve as a valuable resource for pig basic biological discovery, pig genetics, breeding, biotechnology, domestication, veterinary and human biomedicine research.



## 2. Analysis pipelines

This repository contains analysis pipelines used in meta-GWAS, including:

- Genotype imputation and quality control

- Deregressed proofs calculation and individual GWAS analyses
- Meta-level quality control for summary statistics of individual GWASs
- Meta-GWAS analysis and QTL defining
- Calculate the posterior probability that lead SNP has an effect in each study
- QTL validation and genomic prediction on breed
- Functional annotation and enrichment (snpEff, chromatin state, genomic categories, conserved elements, tissue-specific functional regions)
- Genetic parameter estimation based on GWAS summary and heritability enrichment for molecular QTL

- Comparative analysis between pig and human

Pipeline components are public available in this repository and execution scripts are provided in snakemake workflow (.smk).

## 3. Citation

**Integrating large-scale meta-GWAS and PigGTEx resources to decipher the genetic basis of complex traits in pig**

http://dx.doi.org/10.1101/2023.10.09.561393

Zhiting Xu†, Qing Lin†, Xiaodian Cai†, Zhanming Zhong†, Jinyan Teng†, Bingjie Li, Haonan Zeng, Yahui Gao, Zexi Cai, Xiaoqing Wang, Liangyu Shi, Xue Wang, Yi Wang, Zipeng Zhang, Yu Lin, Shuli Liu, Hongwei Yin, Zhonghao Bai, Chen Wei, Jun Zhou, Wenjing Zhang, Xiaoke Zhang, Shaolei Shi, Jun Wu, Shuqi Diao, Yuqiang Liu, Xiangchun Pan, Xueyan Feng, Ruiqi Liu, Zhanqin Su, Chengjie Chang, Qianghui Zhu, Yuwei Wu, The PigGTEx Consortium, Zhongyin Zhou, Lijing Bai, Kui Li, Qishan Wang, Yuchun Pan, Zhong Xu, Xianwen Peng, Shuqi Mei, Delin Mo, Xiaohong Liu, Hao Zhang, Xiaolong Yuan, Yang Liu, George E. Liu, Guosheng Su, Goutam Sahana, Mogens Sandø Lund, Li Ma, Ruidong Xiang, Xia Shen, Pinghua Li, Ruihuang Huang, Maria Ballester, Daniel Crespo-Piazuelo, Marcel Amills, Alex Clop, Peter Karlskov-Mortensen, Merete Fredholm, Guoqing Tang, Mingzhou Li, Xuewei Li, Xiangdong Ding, Jiaqi Li, Yaosheng Chen, Qin Zhang, Yunxiang Zhao *, Fuping Zhao *, Lingzhao Fang *, Zhe Zhang *

bioRxiv, 2023, DOI: [10.1101/2023.10.09.561393](http://dx.doi.org/10.1101/2023.10.09.561393)