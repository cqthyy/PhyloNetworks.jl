using PhyloNetworks
using DataFrames
using StatsModels

n = 4; m = 3
starnet = readTopology(PhyloNetworks.startree_newick(n))
netnames = ["t1","t2","t3","t4"] # = tipLabels(starnet) # was generated to have length n
# individual level
df = DataFrame(
    species = repeat(netnames, inner=m),
    trait1 = [2.8564, missing, 2.8564, 2.8457, 2.8457, 2.8457, 0.4197, 0.4197, 0.4197, 2.2359, 2.2359, 2.2359],
    trait2 = [-2.0935, -2.0935, -2.0935, 0.4955, 0.4955, 0.4955, -2.1977, -2.1977, -2.1977, -2.618, -2.618, -2.618],
    trait3 = [missing, 15.7869, 14.0615, 14.3547, 13.3721, 15.7062, 9.2764, 8.8553, 8.7627, 15.0298, 15.8258, 15.3248]
    # missing replaced 14.0467
)
Y = df[!,:trait3] # nm vector
X = fill(1.0, (n,2)); X[:,2] = df[1:m:n*m,:trait1] # nx2 matrix
#= reduced data: one row per species, extra columns for SD and n of trait 3
gdf = groupby(df, :species)
df_r = combine(gdf, :trait1 => (x -> mean(skipmissing(x))) => :trait1,
                    :trait2 => (x -> mean(skipmissing(x))) => :trait2,
                    :trait3 => (x -> mean(skipmissing(x))) => :trait3,
                    :trait3 => (x ->  std(skipmissing(x))) => :trait3_sd,
                    :trait3 => (x -> sum(.!ismissing.(x))) => :trait3_n)
=#
# as many rows as there are species
df_r = DataFrame(
  species = ["t1","t2","t3","t4"], trait1 = [2.8564,2.8457,.4197,2.2359],
  trait2 = [-2.0935,0.4955,-2.1977,-2.618],
  trait3 = [14.0615, 14.477666666666666, 8.9648, 15.393466666666667],
  #        [14.9242,... if we include ind 2
  trait3_sd = [0,1.1718985891848046,.2737966581242371,.4024181076111425],
  # [1.2200420402592675,... if we include ind 2
  trait3_n = [1,3,3,3]) # [2,3,3,3]

m1 = phylolm(@formula(trait3 ~ trait1), df, starnet; reml=true,
      tipnames=:species, withinspecies_var=true)
m2 = phylolm(@formula(trait3 ~ trait1), df, starnet; model="lambda", reml=true,
      tipnames=:species, withinspecies_var=true, fixedValue=1)