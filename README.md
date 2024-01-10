# hybrid2simulation

## Overview

This repository contains all the necessary R code to run a simulation study that compares the various study design methods available for hybrid type 2 CRTs. It uses the 'hybrid2power' R package, which can be found [here](https://github.com/melodyaowen/hybrid2power/tree/main). 

## Aims

The specific aims of the simulation study are to evaluae the impacts on the estimands of
1. constricting the number of clusters, $K$, to varying values (particularly small values);
2. of varying the cluster sizes, $m$, while still assuming equal cluster sizes among all clusters;
3. of varying the effect size ($\beta_1$ and $\beta_2$), endpoint specific ICCs ($\rho_0^{(1)}$ and $\rho_0^{(2)}$), and total outcome variances ($\sigma_1^2$ and $\sigma_2^2$);
4. of varying the intra-subject between-endpoint ICC ($\rho_1^{(1,2)}$) and the inter-subject between-endpoint ICC ($\rho_2^{(1,2)}$);
5. of changing the outcome types, i.e. either two continuous endpoints or two binary endpoints;
6. of varying the treatment allocation ratio, $r$

## Data-Generating Mechanisms

Two key cases will be examined resulting in two data-generating mechanism schemes - the case of two co-primary continuous outcomes, and two co-primary binary outcomes.

## Estimands

The key estimands that will be examined are the following:
1. Statistical power, $\pi$
2. Required number of clusters, $K$
3. Cluster size, $m$
4. Type I error
5. Minimum detectable effect plane for $\beta_1$ and $\beta_2$
