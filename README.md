demixed Principal Component Analysis (dPCA)
===========================================

dPCA is a linear dimensionality reduction technique that automatically discovers and highlights the essential features of complex population activities. The population activity is decomposed into a few demixed components that capture most of the variance in the data and that highlight the dynamic tuning of the population to various task parameters, such as stimuli, decisions, rewards, etc.

> D Kobak<sup>+</sup>, W Brendel<sup>+</sup>, C Constantinidis, CE Feierstein,
A Kepecs, ZF Mainen, X-L Qi, R Romo, N Uchida, CK Machens<br>
> **Demixed principal component analysis of neural population data**<br>
> eLife 2016, https://elifesciences.org/content/5/e10989<br>
> (arXiv link: http://arxiv.org/abs/1410.6031)

This repository provides preprocessing and some of the analysis scripts for the dPCA paper. They allow to reproduce Figures 3 to 6 of the paper (the ones yielding dPCA components of all four datasets), starting with the raw data available at:

1. https://crcns.org/data-sets/pfc/pfc-4
2. https://crcns.org/data-sets/pfc/pfc-3
3. https://crcns.org/data-sets/ofc/ofc-1
4. https://crcns.org/data-sets/ofc/ofc-2

The dPCA toolkit is located at https://github.com/machenslab/dPCA.

-----------

**Note (May 2016):** The scripts are currently provided as is, but will be adapted during summer 2016 to be more readable.

**Update (Jul 15, 2016):** The scripts processing the Romo dataset are ready.
