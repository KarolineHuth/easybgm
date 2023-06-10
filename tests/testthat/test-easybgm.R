suppressPackageStartupMessages({
  library(bgms)
  library(vdiffr)
})

##--------------------------------
## Fitting with bgms
##--------------------------------

set.seed(123)
data <- na.omit(Wenchuan)
res_bgms <- easybgm(data[1:100, 1:5], type = "ordinal",
                    package = "bgms", save = T)
test_that("easybgm works for bgms", {
  expect_snapshot(summary(res_bgms))
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
poststruc_bgms <- plot_posteriorstructure(res_bgms)
vdiffr::expect_doppelganger("posterior structure plot bgms", poststruc_bgms)

# 4. Posterior complexity plot
postcompl_bgms <- plot_posteriorcomplexity(res_bgms)
vdiffr::expect_doppelganger("posterior complexity plot bgms", postcompl_bgms)

# 5. structure plot
struc_bgms <- plot_structure(res_bgms)
vdiffr::expect_doppelganger("structure plot bgms", struc_bgms)

# 6. HDI plot
HDI_bgms <-plot_parameterHDI(res_bgms)
vdiffr::expect_doppelganger("HDI plot bgms", HDI_bgms)

##--------------------------------
## Fitting with BDgraph
##--------------------------------
set.seed(123)
res_bdgraph <- easybgm(data[1:100, 1:5], type = "continuous",
                    package = "BDgraph", save = T)
test_that("easybgm works for bdgraph", {
  expect_snapshot(summary(res_bdgraph))
})

##--------------------------------
## Plotting with bgms
##--------------------------------


##--------------------------------
## Fitting with BGGM
##--------------------------------
# DOES NOT WORK, output keeps changing slightly despite set.seed
# set.seed(123)
# res_bggm <- easybgm(data[1:300, 1:5], type = "continuous",
#                        package = "BGGM")
# test_that("easybgm works for bggm", {
#   expect_snapshot(summary(res_bggm))
# })



