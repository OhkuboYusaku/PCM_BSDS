# "BSDSモデルによる方向性淘汰圧の推定"
author: "Ph.D. Ohkubo Yusaku (ROIS-DS & Institute of Statistical Mathematics) {y-ohkubo[--]ism.ac.jp}"<br>
date: "May.25 2022"

枝特異的方向性淘汰モデル(branch-specific directional selection; 以下、BSDS)は連続形質に用いる種間系統比較法(Phylogenetic Comparative Method)の一種で、系統樹における一部の枝で方向性淘汰を経た生物の形質を分析するために開発されました。このページでは、BSDSモデルを用いる３つの推定方法について紹介します。


The branch-specific directional selection model (BSDS) is a model of Phylogenetic Comparative Methods  (PCM) for a continuous trait. The authors developed this model to analyze phenotypic traits that experienced directional selection in a specific edge of a whole phylogenetic tree. Here on this page, I introduce three different applications of the BSDS model. 

__各種の平均形質値を使って直接多変量正規分布を当てはめる場合は[こちら](https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_MLE)を参照してください。最尤法を用いて祖先形質(MRCA)、進化率(ev)、選択圧の強さ(k)などマクロ進化に関わるパラメータを推定する方法(Rコード)が紹介されています。__ この方法はもっとも計算速度が速く大きな種数(例えば1000種)のデータや、個体毎の形質が手に入らない場合にもにも適用することができます。しかし各種の平均形質値が既知として扱われるために種内分散や平均値の推定誤差は考慮されません。また、種数が少ない(例えば~50種)場合にはevやkの信頼区間が甘めに見積もられる場合があります。

__A simple multivariate-normal distribution approach is described [here](https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_MLE), which fits mean trait values. This page shows how to estimate the most common recent ancestor (MRCA), evolutionary rate (ev), and the strength of the directional selection using the maximum likelihood estimator, including the R-code.__ This approach is the most computationally efficient methodology and scales to a larger dataset (e.g. N_sp=1000). It applies even when we do not have individual trait data. But, it dismisses estimation error of the mean trait values and intra-species variations. Further, the confidence interval of ev, k could be too optimistic when the number of species is small (e.g. <50).

階層ベイズモデルによって、各種の平均形質とPCMのパラメータ(MRCA, ev, k)を同時に推定する方法(RおよびStanコード)はこちらを参照してください。
この方法は、平均形質の推定誤差や種内分散を補正しながら(MRCA, ev, k)を推定できますが、計算時間は大きくなります(10サンプルずつ50種の場合でおよそ1-2分)。

A hierarchical Bayesian modeling approach is described here, including the R and Stancode. Simultaneous estimation of both a mean trait value within each species and the latent parameter behind the evolutionary history is possible. 
This approach adjusts estimation error and intra-species variation when we estimate (MRCA, ev, k), but computationally more demanding (a few minutes for 10 samples * 50 species dataset).

__BSDSモデルをランダム効果回帰モデル(BSDS-LMM)としてあつかい、回帰係数などを推定する方法(RとStanコード)は[こちら](https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_LMM)を参照してください__。
複数種から得られたデータを分析する際に、進化的な近縁度によって生じた”みせかけの相関”をBSDSで補正しながら説明変数Xの回帰係数を推定することができます。

__The linear mixed model, which assumes BSDS for the covariance structure of random-effect(BSDS-LMM) is described [here](https://github.com/OhkuboYusaku/PCM_BSDS/tree/main/example/BSDS_LMM)__. 
When a dataset contains multiple species, this model corrects "spurious correlation", induced by the shared evolutionary history among species,  when we estimate a regression coefficient of an explanatory variable X.
