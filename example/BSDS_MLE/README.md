
## Introduction
枝特異的方向性淘汰モデル(branch-specific directional selection; 以下、BSDS)は種間系統比較法(phylogenetic comparative method)の一種で、系統樹における一部の枝で方向性淘汰を経た生物の形質を分析するために開発されました。このドキュメントでは、仮想の形質データを題材にRで実際にBSDSモデルを使った分析を実行する方法について紹介します。


The branch-specific directional selection model (BSDS) is one of the statistical methodologies called the phylogenetic comparative methods (PCM). It aims to analyze cases where focal phenotypic traits experienced “directional selection” at only parts of a phylogenetic tree. In addition to the most recent common ancestor (MRCA), evolution rate, this model estimates the strength of the directional selection at a given edge as well as its uncertainty. This document introduces a practical guide to conducting the maximum likelihood inferences using multivariate normal distributions.



__種内変動を考慮した階層モデルにより推定したい場合は、以下を参照してください。__

__This page assumes that the true trait mean is known (e.g. observation error does not exist or is almost negligible). To take within-species variation into account, see the hierarchical modeling approach below.__

https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_LMM

## 下準備/preparations
### 必要パッケージのインストール/Installing necessarily packages
以下のコードでは、{ape},{mvnfast}の2つのパッケージに依存しています。事前にインストールし、読み込んでおきます。

The following example program depends on two R packages, {ape} and {mvnfast}. Make sure to be installed, and import in advance.

```r
install.packages(c("ape", "mvnfast"))
```


```r
library(ape)
library(mvnfast)
```


### 必要な関数の準備/ definition of the necessarily function
BSDSモデルで最尤推定するには、独自の尤度関数を定義しておく必要があります。

To obtain the maximum likelihood estimate, define the log-likelihood function in advance.

```r
 loglik<- function(phylo, theta, y, DS_edge, MRCA_ij){
    MRCA<- theta[1]
    ev<- (theta[2])
    k<- numeric(length(phylo$edge.length)) + 1
    for(i_DS in 1:length(DS_edge)){
      k[DS_edge[i_DS]]<- exp(theta[3])
    }
    
    #define data
    tree_obj<- as.matrix(phylo[["edge"]])
    simulated_mean_trait<- numeric(length(phylo$edge.length))
    simulated_var_trait<- numeric(length(phylo$edge.length))
    len_phylo<- length(phylo$edge.length)
    ###### NOTE ######
    # simulated_mean_trait[i] and simulated_var_trait[i] contains simulated trait of sp i.
    
    ## conduct evolution simulation
    simulated_mean_trait[tree_obj[1,1]]<- MRCA
    simulated_var_trait[tree_obj[1,1]]<- 0
    
    for(i in 1:len_phylo){
      branch_len=phylo$edge.length[i]
      simulated_mean_trait[tree_obj[i,2]]<- simulated_mean_trait[tree_obj[i,1]] + (branch_len*ev*(k[i]^2-1))/k[i]
      simulated_var_trait[tree_obj[i,2]]<- simulated_var_trait[tree_obj[i,1]] + (2*branch_len*ev*(k[i]^2+1))/k[i]
    }
    
    # output vcov matrix ot tip sps.
    N_tip<- len_phylo - phylo$Nnode +1
    vcv_BSDE<- matrix(0, N_tip, N_tip)
    
    for(i in 1: N_tip){
      for(j in i: N_tip){
        if(i != j){
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[MRCA_ij[i,j]]
        }else{
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[i]
        }
      }
    }

    loglik<- try(dmvn(y, simulated_mean_trait[1:N_sp],  Matrix::nearPD(vcv_BSDE)$mat, log=T))
    if(is.numeric(loglik)==F){return(-10000)}
    
    return(loglik)
  }
```


## 実行例/running example
### データの読み込み/ import data
まず、題材となるデータを読み込み、構造を確認します。

First, let us import data and check its structure.
```r
data<- read.csv("BSDS_MLE.csv")
summary(data)
```

```
##        Y              sp_ID     
##  Min.   : 91.14   Min.   : 1.0  
##  1st Qu.:100.10   1st Qu.: 3.0  
##  Median :104.75   Median : 5.5  
##  Mean   :104.43   Mean   : 5.5  
##  3rd Qu.:107.02   3rd Qu.: 8.0  
##  Max.   :119.54   Max.   :10.0
```

```r
Y<- data$Y
sp_ID<- (data$sp_ID)
```

Yに各種(i=1,2,...N_sp)の形質値を格納しています。

Y contains trait data of each individuals (i=1, 2, …N_sample),and sp_ID is an indicator of species ID (1,2, …N_sp).

次に系統樹を読み込み加工しておきます。ここでは、{ape}パッケージの関数を用いてNewick形式で記録された系統樹を読み込みます。また、二つの種i, jの共通祖先を定義しておきます。

Next, let us import a phylogenetic tree. Here in this example, we use a function from {ape} package and import a Newick formatted file. The common ancestor of the two species i,j is also defined. 
```r
phylo<- read.tree("BSDS_LMM_tree")
plot(phylo)
axisPhylo()
 MRCA_ij<- matrix(0, N_tip, N_tip) ## i,j elements correspond to the location of their MRCA in the tree
 
  for(i in 1:N_tip){
    for(j in i:N_tip){
      MRCA_ij[i,j]<- MRCA_ij[j,i]<- getMRCA(phylo, tip=c(i,j))
    }
  }
```

![](BSDS_MLE_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
phylo$edge # tree構造:[,1]の親種から[,2]の子孫種へエッジが伸びている
# tree structure: Edges are connected from a parent species [,1] to a descendent species [,2]

```

```
##       [,1] [,2]
##  [1,]   11   12
##  [2,]   12   13
##  [3,]   13    1
##  [4,]   13    2
##  [5,]   12   14
##  [6,]   14   15
##  [7,]   15   16
##  [8,]   16   17
##  [9,]   17    3
## [10,]   17    4
## [11,]   16    5
## [12,]   15    6
## [13,]   14    7
## [14,]   11   18
## [15,]   18    8
## [16,]   18   19
## [17,]   19    9
## [18,]   19   10
```

また、方向性淘汰が生じたと思われる箇所を指定します。例えば、「原生種t10が、t1, t6, t9の共通祖先(sp19)と分岐した直後から方向性淘汰を受けるようになった」と仮定します。この場合、phylo$edgeにおいて「sp19からsp10」に伸びる枝、すなわち[18,]が該当する箇所になります。


```r
DS_edge<- 18
```


最後に、これらのデータから種ごとの平均値を求めておきます。

Finally, obtain the mean trait for each species. 

```r
  model<- lm(Y~as.factor(sp_ID)+0)
```


### 推定の実行/conduct inferences
R標準のoptim関数で尤度関数の最大化を行います。

Maximize the log-likelihood function using optim function.

```r
  BSDS_MLE<- optim(par=c(100, 10, 5), fn=loglik, gr = NULL,
                       phylo=phylo, y=model$coefficients, DS_edge=DS_edge,
                       method = "BFGS", MRCA_ij=MRCA_ij,
                       control=list(fnscale=-1, trace=0, maxit=1000000), hessian=T)
```

parに、MRCA, ev, kの適当な初期値を与えてあります。うまくいかない場合も初期値をかえるか、事前にmethod = "SANN" (シミュレーテッドアニーリング法のこと)を指定して”予熱”しておいてからBFGSに渡すと正常に動作する場合があります。

Here, par is an initial value of MRCA, ev rate, and k, the strength of directional selection, to be evaluated. If failed to optimize, take different initial values or "preheat" the log-likelihood using method =" SANN", prior to BFGS method.

## 結果の出力/output
最尤推定値、その標準誤差、Z値などを抽出できます。なお、SEの評価には素朴に尤度関数のヘッセ行列を使っています(即ち、最尤推定量に関して通常の漸近理論を仮定)。しかし、データの種数が少ない場合(例えば10種)には見積もりがあまくなる場合があるようです。ベイズ推定を用いることで見積もりが改善する場合があります。

You can extract the maximum likelihood estimate, its standard error, Z-value, etc. Here, the Hessian matrix is used to evaluate SE (i.e. the standard asymptotic theory of MLE is assumed). It tends to be, however, too optimistic when the number of species is too small (e.g. N_sp ＝10). Bayes estimate could be help in this case.

```r
MLE<- BSDS_MLE$par[1:3] #MLE
SE<- sqrt(diag(solve(-BSDS_MLE$hessian))) # SE
Z<- MLE/SE　# Z-value

CI_lo<- MLE - 1.96*SE
CI_up<- MLE + 1.96*SE
```


## References
Ohkubo, Kutsukake, Koizumi (under review) "Evaluating a strength of directional selection using a novel branch-specific directional selection (BSDS) model of phylogenetic comparative method"
