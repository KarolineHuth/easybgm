suppressPackageStartupMessages(library(bgms))
data <- na.omit(Wenchuan)

##--------------------------------
## Fitting with bgms
##--------------------------------
set.seed(123)
res_bgms <- easybgm(data[1:100, 1:5], type = "ordinal",
                    package = "bgms", save = T, centrality = T, iter = 1000)
test_that("easybgm works for bgms", {
  testthat::expect_snapshot(summary(res_bgms))
})

##--------------------------------
## Plotting with bgms
##--------------------------------

# 1. Network plot
network_bgms <- plot_network(res_bgms)
vdiffr::expect_doppelganger("network plot bgms", network_bgms)

# 2. Evidence plot
evidence_bgms <- plot_network(res_bgms)
vdiffr::expect_doppelganger("evidence plot bgms", evidence_bgms)

# 3. Posterior structure plot
poststruc_bgms <- plot_structure_probabilities(res_bgms)
vdiffr::expect_doppelganger("posterior structure plot bgms", poststruc_bgms)

# 4. Posterior complexity plot
postcompl_bgms <- plot_complexity_probabilities(res_bgms)
vdiffr::expect_doppelganger("posterior complexity plot bgms", postcompl_bgms)

# 5. structure plot
struc_bgms <- plot_structure(res_bgms)
vdiffr::expect_doppelganger("structure plot bgms", struc_bgms)

# 6. HDI plot
HDI_bgms <-plot_parameterHDI(res_bgms)
vdiffr::expect_doppelganger("HDI plot bgms", HDI_bgms)

# 7. centrality plot
centrality_bgms <-plot_centrality(res_bgms)
vdiffr::expect_doppelganger("centrality plot bgms", centrality_bgms)

##--------------------------------
## Fitting with BDgraph
##--------------------------------
set.seed(123)
res_bdgraph <- suppressWarnings(easybgm(data[1:100, 1:5], type = "continuous",
                    package = "BDgraph", save = T, centrality = T, iter = 1000))
test_that("easybgm works for bdgraph", {
  testthat::expect_snapshot(summary(res_bdgraph))
})

##--------------------------------
## Plotting with BDgraph
##--------------------------------

#1. Network plot
network_bdgraph <- plot_network(res_bdgraph)
vdiffr::expect_doppelganger("network plot Bdgraph", network_bdgraph)

# 2. Evidence plot
evidence_bdgraph <- plot_network(res_bdgraph)
vdiffr::expect_doppelganger("evidence plot Bdgraph", evidence_bdgraph)

# 3. Posterior structure plot
poststruc_bdgraph <- plot_structure_probabilities(res_bdgraph)
vdiffr::expect_doppelganger("posterior structure plot Bdgraph", poststruc_bdgraph)

# 4. Posterior complexity plot
postcompl_bdgraph <- plot_complexity_probabilities(res_bdgraph)
vdiffr::expect_doppelganger("posterior complexity plot Bdgraph", postcompl_bdgraph)

# 5. structure plot
struc_bdgraph <- plot_structure(res_bdgraph)
vdiffr::expect_doppelganger("structure plot Bdgraph", struc_bdgraph)

# 6. HDI plot
HDI_bdgraph <- suppressWarnings(plot_parameterHDI(res_bdgraph))
vdiffr::expect_doppelganger("HDI plot Bdgraph", HDI_bdgraph)

# 7. centrality plot
centrality_bdgraph <-plot_centrality(res_bdgraph)
vdiffr::expect_doppelganger("centrality plot Bdgraph", centrality_bdgraph)

##--------------------------------
## Fitting with BGGM
##--------------------------------
# DOES NOT WORK, output keeps changing slightly despite set.seed
# set.seed(123)
# res_bggm <- easybgm(data[1:300, 1:5], type = "continuous",
#                        package = "BGGM")
# test_that("easybgm works for bggm", {
#   vdiffr::expect_snapshot(summary(res_bggm))
# })


##--------------------------------
## Fitting with bgms
##--------------------------------


set.seed(123)
data <- na.omit(Wenchuan)
fit_bgms <- bgm(data[1:100, 1:5], iter = 1000)
network_bgmfit <- plot_network(fit_bgms)
vdiffr::expect_doppelganger("network plot using bgm to fit", network_bgmfit)


##--------------------------------
## Sparse vs dense test
##--------------------------------
set.seed(123)
data <- na.omit(Wenchuan)
sparse_dense <- sparse_or_dense(data[1:100, 1:5], type = "ordinal", iter = 1000)
test_that("easybgm works for sparse vs dense", {
  testthat::expect_snapshot(sparse_dense)
})
