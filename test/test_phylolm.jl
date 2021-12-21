using PhyloNetworks
using DataFrames
using StatsModels

# "phylolm: binary tree, withinspecies var"
net = readTopology("(((t4:0.2537636499,t5:0.2537636499):0.2103870459,t3:0.4641506959):0.5358493041,(t1:0.4807642475,t2:0.4807642475):0.5192357525);")
df = DataFrame(
      y=[11.399108,3.216645,13.648011,4.851454,-4.922803],
      y_sd=[0.2463003,0.3236629,0.2458547,0.4866844,0.3434582],
      y_n=fill(5,5),
      x1=[0.7894387,0.3285286,2.1495131,-1.4935665,-2.3199935],
      x2=[2.8184268,0.4740687,2.5801004,1.9967379,-0.4121874],
      species=net.names
)

n = 5
starnet = readTopology(PhyloNetworks.startree_newick(n))

# for 0 and 1 - can compare results to ordinary BM model
# lambda = 0, compare to star tree
m5 = phylolm(@formula(y ~ x1), df, starnet; reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)
#=
StatsModels.TableRegressionModel{PhyloNetworkLinearModel, Matrix{Float64}}

Formula: y ~ 1 + x1

Model: Brownian motion

Parameter Estimates, using REML:
phylogenetic variance rate: 15.5623
within-species variance: 0.116138

Coefficients:
──────────────────────────────────────────────────────────────────────
               Coef.  Std. Error     t  Pr(>|t|)  Lower 95%  Upper 95%
──────────────────────────────────────────────────────────────────────
(Intercept)  6.03325     1.76961  3.41    0.0422   0.401575   11.6649
x1           3.61458     1.09896  3.29    0.0461   0.117202    7.11195
──────────────────────────────────────────────────────────────────────
Log Likelihood: -21.3324812448
AIC: 50.6649624896
=#

m4 = phylolm(@formula(y ~ x1), df, net; model="lambda", reml=true,
      tipnames=:species, withinspecies_var=true, fixedValue=0.0, y_mean_std=true)
#=
StatsModels.TableRegressionModel{PhyloNetworkLinearModel, Matrix{Float64}}

Formula: y ~ 1 + x1

Model: Pagel's lambda

Parameter Estimates, using REML:
phylogenetic variance rate: 15.5623
Lambda: 0
within-species variance: 0.116138

Coefficients:
──────────────────────────────────────────────────────────────────────
               Coef.  Std. Error     t  Pr(>|t|)  Lower 95%  Upper 95%
──────────────────────────────────────────────────────────────────────
(Intercept)  6.03325     1.76961  3.41    0.0422   0.401575   11.6649
x1           3.61458     1.09896  3.29    0.0461   0.117202    7.11195
──────────────────────────────────────────────────────────────────────
Log Likelihood: -21.3324812448
AIC: 52.6649624896
=#


# ordinary BM
m3 = phylolm(@formula(y ~ x2), df, net; reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)
#=
StatsModels.TableRegressionModel{PhyloNetworkLinearModel, Matrix{Float64}}

Formula: y ~ 1 + x2

Model: Brownian motion

Parameter Estimates, using REML:
phylogenetic variance rate: 9.30581
within-species variance: 0.116126

Coefficients:
─────────────────────────────────────────────────────────────────────────
                 Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────────
(Intercept)  -0.803512    2.11612   -0.38    0.7294   -7.53794    5.93092
x2            4.1111      0.702594   5.85    0.0100    1.87513    6.34707
─────────────────────────────────────────────────────────────────────────
Log Likelihood: -19.6834453226
AIC: 47.3668906453
=#

m2 = phylolm(@formula(y ~ x2), df, net; model="lambda", reml=true,
      tipnames=:species, withinspecies_var=true, fixedValue=1.0, y_mean_std=true)
#=
StatsModels.TableRegressionModel{PhyloNetworkLinearModel, Matrix{Float64}}

Formula: y ~ 1 + x2

Model: Pagel's lambda

Parameter Estimates, using REML:
phylogenetic variance rate: 9.30581
Lambda: 1
within-species variance: 0.116126

Coefficients:
─────────────────────────────────────────────────────────────────────────
                 Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────────
(Intercept)  -0.803512    2.11612   -0.38    0.7294   -7.53794    5.93092
x2            4.1111      0.702594   5.85    0.0100    1.87513    6.34707
─────────────────────────────────────────────────────────────────────────
Log Likelihood: -19.6834453226
AIC: 49.3668906453
=#

# lambda = 0.5
m1 = phylolm(@formula(y ~ x1), df, net; model="lambda", reml=true,
      tipnames=:species, withinspecies_var=true, fixedValue=0.5, y_mean_std=true)
#=
StatsModels.TableRegressionModel{PhyloNetworkLinearModel, Matrix{Float64}}

Formula: y ~ 1 + x1

Model: Pagel's lambda

Parameter Estimates, using REML:
phylogenetic variance rate: 21.8812
Lambda: 0.5
within-species variance: 0.11614

Coefficients:
──────────────────────────────────────────────────────────────────────
               Coef.  Std. Error     t  Pr(>|t|)  Lower 95%  Upper 95%
──────────────────────────────────────────────────────────────────────
(Intercept)  6.2069      2.5429   2.44    0.0924   -1.88573   14.2995
x1           3.93469     1.42465  2.76    0.0700   -0.59919    8.46858
──────────────────────────────────────────────────────────────────────
Log Likelihood: -21.4050221575
AIC: 52.810044315
=#

# optimize lambda case (where fixedValue is not provided)
m6 = phylolm(@formula(y ~ x1 + x2), df, net; model="lambda", reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)

m7 = phylolm(@formula(y ~ x1 + x2), df, net; reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)