# Estimating Directional Selection using the BSDS Model

**Author:** Dr. Yusaku Ohkubo (ROIS-DS & Institute of Statistical Mathematics)\
**Contact:** y-ohkubo[--]okayama-u.ac.jp\
**Updated:** July 9, 2025

---

## 🔍 Introduction

The Branch-Specific Directional Selection (BSDS) model is a type of phylogenetic comparative method designed to detect directional selection acting on traits along specific branches of a phylogenetic tree.

This document demonstrates how to apply the BSDS model to simulated trait data using R and Stan.

> **Note:** For the BSDS model implemented as a random effect in regression models (BSDS-LMM), please refer to:\
> [https://github.com/OhkuboYusaku/PCM\_BSDS/tree/main/example/BSDS\_LMM](https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_LMM)

---

##  Preparation

### Required Packages

Install and load the following R packages:

```r
install.packages(c("ape", "rstan", "dummies"))

library(ape)
library(rstan)
```

If you encounter issues installing Stan, refer to:

- Japanese guide: [https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started-(Japanese)](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started-\(Japanese\))
- English guide: [https://mc-stan.org/users/interfaces/rstan](https://mc-stan.org/users/interfaces/rstan)

---

### Data Preparation Function

To run the BSDS model in Stan, trait values, species IDs, and the phylogenetic tree (topology and branch lengths) must be passed as a list. The following function prepares this data:

```r
BSDS2stan_data <- function(phylo, y, Z, D_edge) {
  len_phylo <- length(phylo$edge.length)
  N_tip <- len_phylo - phylo$Nnode + 1
  branch_len <- phylo$edge.length
  tree_obj <- as.matrix(phylo$edge)
  MRCA_ij <- matrix(0, N_tip, N_tip)

  for (i in 1:N_tip) {
    for (j in i:N_tip) {
      MRCA_ij[i, j] <- MRCA_ij[j, i] <- getMRCA(phylo, tip = c(i, j))
    }
  }

  list(
    N = length(y), N_sp = N_tip, len_phylo = len_phylo,
    branch_len = branch_len, tree_obj = tree_obj, MRCA_ij = MRCA_ij,
    y = y, Z = dummies::dummy(Z), DS_edge = D_edge
  )
}
```

---

##  Running the Model

### 1. Load Trait Data

```r
data <- read.csv("BSDS_sample.csv")
summary(data)

Y <- data$Y
sp_ID <- data$sp_ID
```

### 2. Load and Inspect the Phylogenetic Tree

```r
phylo <- read.tree("BSDS_tree")
plot(phylo)
axisPhylo()
phylo$edge  # Display tree structure
```

### 3. Specify the Branch Under Directional Selection

For example, if directional selection begins right after the divergence from a common ancestor between species t1 and t5:

```r
DS_edge <- 2  # Refers to the 2nd row of phylo$edge
```

### 4. Prepare Data for Stan

```r
dat <- BSDS2stan_data(phylo, Y, sp_ID, DS_edge)
```

---

##  Stan Settings and Sampling

```r
scr <- "stan_BSDS.stan"
war <- 5000
ite <- 25000
cha <- 2
par <- c("MRCA", "ev", "sel", "log_likelihood")

options(mc.cores = parallel::detectCores())
```

```r
fit_BSDS <- stan(
  file = scr, model_name = scr, data = dat, pars = par,
  chains = cha, warmup = war, iter = ite, thin = 10
)
```

> **Note:** If you encounter warnings such as "divergent transitions" or "maximum treedepth exceeded," consult:\
> [https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup](https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup)

---

## Inspecting Results

Ensure convergence by checking trace plots and ˆR values:

```r
print(fit_BSDS)
traceplot(fit_BSDS)
```

If available, use the **{shinystan}** package for interactive diagnostics and visualization.

---

## References

- Ohkubo, Y., Kutsukake, N., & Koizumi, I. (2023).\
  *A novel phylogenetic comparative method for evaluating the strength of branch-specific directional selection.* Evolution, 77(1), 63–82.

- Ohkubo, Y., Kutsukake, N., & Koizumi, I. (2021).\
  "Estimating Branch-Specific Directional Selection," presented at the 68th Annual Meeting of the Ecological Society of Japan.\
  [https://esj.ne.jp/meeting/abst/68/D01-12.html](https://esj.ne.jp/meeting/abst/68/D01-12.html)

