---
title: "単純な多変量正規BSDSモデルを使った最尤推定"
author: "Ph.D. Ohkubo Yusaku (ROIS-DS & Institute of Statistical Mathematics) {y-ohkubo[--]ism.ac.jp}"
date: "2/Dec./2021"
output: 
  html_document: 
    keep_md: yes
    self_contained: TRUE
 
---



## Introduction
枝特異的方向性淘汰モデル(branch-specific directional selection; 以下、BSDS)は種間系統比較法(phylogenetic comparative method)の一種で、系統樹における一部の枝で方向性淘汰を経た生物の形質を分析するために開発されました。
このドキュメントでは、仮想の形質データを題材にRで実際にBSDSモデルを使った分析を実行する方法について紹介します。


__種内変動を考慮した階層モデルにより推定したい場合は、以下を参照してください。__


## 下準備
### 必要パッケージのインストール
以下のコードでは、{ape},{mvnfast}{dummies}の3つのパッケージに依存しています。事前にインストールし、読み込んでおきます。

```r
install.packages(c("ape", "mvnfast"))
```


```r
library(ape)
library(mvnfast)
```


### 必要な関数の準備
BSDSモデルで最尤推定するには、独自の尤度関数を定義しておく必要があります。

```r
# define data-arrangement function
  loglik<- function(phylo, theta, y, DS_edge){
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
          MRCA_ij<- getMRCA(phylo, tip=c(i,j))
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[MRCA_ij]
        }else{
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[i]
        }
      }
    }
    
    if((sum(is.na(vcv_BSDE))>0)||(min(eigen(vcv_BSDE)$values)<=0)){return(-Inf)}
    loglik<- try(dmvn(y, simulated_mean_trait[1:N_tip],  Matrix::nearPD(vcv_BSDE)$mat, log=T))
    if(is.numeric(loglik)==F){return(-10000)}
    
    return(loglik)
  }
```


## 実行例
### データの読み込み
まず、題材となるデータを読み込み、構造を確認します。


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

Yに各個体(i=1,2,...N_sample)の形質値、sp_IDに各個体の種ID(1,2,...N_sp)を格納しています。


次に系統樹を読み込み加工しておきます。ここでは、{ape}パッケージの関数を用いてNewick形式で記録された系統樹を読み込みます。

```r
phylo<- read.tree("BSDS_LMM_tree")
plot(phylo)
axisPhylo()
```

![](BSDS_MLE_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
phylo$edge # tree構造:[,1]の親種から[,2]の子孫種へエッジが伸びている
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


```r
  model<- lm(Y~as.factor(sp_ID)+0)
```


### 推定の実行
R標準のoptim関数で尤度関数の最大化を行います。

```r
  BSDS_MLE<- optim(par=c(100, 10, 5), fn=loglik, gr = NULL,
                       phylo=phylo, y=model$coefficients, DS_edge=DS_edge,
                       method = "BFGS",
                       control=list(fnscale=-1, trace=0, maxit=1000000), hessian=T)
```

parに、MRCA, ev, kの適当な初期値を与えてあります。うまくいかない場合も初期値をかえるか、事前にmethod = "SANN"を指定して”予熱”しておいてからBFGSに渡すと正常に動作する場合があります。

## 結果の出力
最尤推定値、その標準誤差、Zなどを抽出できます。なお、SEの評価には素朴に尤度関数のヘッセ行列を使っていますがデータの種数が少ない場合(例えば10種)には見積もりがあまくなる場合があるようです。

```r
MLE<- BSDS_MLE$par[1:3] #MLE
SE<- sqrt(diag(solve(-BSDS_MLE$hessian))) # SE
Z<- MLE/SE　# Z-value

CI_lo<- MLE - 1.96*SE
CI_up<- MLE + 1.96*SE
```


## References
Ohkubo, Kutsukake, Koizumi (under review) "Evaluating a strength of directional selection using a novel branch-specific directional selection (BSDS) model of phylogenetic comparative method"
