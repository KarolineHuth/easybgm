#### save parameter checks  ####
### Check centrality computation; see centrality argument requires save = T
### test type = blume-capel for bgms
### prior for bgms SBM
### how do i vary the versions of bgms with easybgm 
devtools::load_all()
test_that("easybgm returns expected structure across valid type–package combos", {
  set.seed(123)
  
  # Subsample small data to stay fast on CRAN
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p <- ncol(dat)
  
  # Test only core combinations
  combos <- list(
    list(type = "continuous", pkg = "BGGM"),
    list(type = "mixed",      pkg = "BDgraph"),
    list(type = "binary",     pkg = "bgms")
  )
  
  for (cmb in combos) {
    t <- cmb$type
    pkg <- cmb$pkg
    
    not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL
    
    suppressMessages({
      res <- easybgm(
        data       = dat,
        type       = t,
        package    = pkg,
        iter       = 10,          # tiny for speed
        save       = FALSE,
        centrality = FALSE,
        progress   = FALSE,
        not_cont   = not_cont
      )
    })
    
    # --- class check ---
    expect_true(inherits(res, c("easybgm")))
    expect_true(any(grepl("package_", class(res))))  # backend tag present
    
    # --- field presence check ---
    expect_true(all(c("parameters", "inc_probs", "inc_BF", "structure", "model") %in% names(res)))
    
    # --- dimensions check ---
    expect_equal(dim(res$parameters), c(p, p))
    expect_equal(dim(res$inc_probs),  c(p, p))
    expect_equal(dim(res$inc_BF),     c(p, p))
    expect_equal(dim(res$structure),  c(p, p))
    
    # --- sanity check ---
    expect_false(all(is.na(res$parameters)))
    expect_false(all(is.na(res$inc_probs))) 
    
    ### if save = TRUE also check for `samples_posterior`with dimension k (= p*(p-1)/2) x iter and `centrality` with dimension p x iter 
    ### if is SBM check for list element sbm with four elements
  }
})

test_that("plotting functions work across valid type–package combos", {
  set.seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p   <- ncol(dat)
  
  combos <- list(
    list(type = "continuous", pkg = "BGGM"),
    list(type = "mixed",      pkg = "BDgraph"),
    list(type = "binary",     pkg = "bgms")
  )
  
  for (cmb in combos) {
    t   <- cmb$type
    pkg <- cmb$pkg
    not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL
    
    suppressMessages({
      res <- easybgm(
        data       = dat,
        type       = t,
        package    = pkg,
        iter       = 10,
        save       = TRUE,
        centrality = TRUE,
        progress   = FALSE,
        not_cont   = not_cont
      )
    })
    
    # --- edge evidence ---
    g1 <- invisible(plot_edgeevidence(res))
    expect_true(inherits(g1, c("ggplot", "qgraph")))
    
    # --- network ---
    g2 <- invisible(plot_network(res))
    expect_true(inherits(g2, c("ggplot", "qgraph")))
    
    # --- structure plots (skip for BGGM) ---
    if (pkg != "BGGM") {
      g3 <- invisible(plot_structure_probabilities(res))
      expect_s3_class(g3, "ggplot")
      
      g4 <- invisible(plot_complexity_probabilities(res))
      expect_s3_class(g4, "ggplot")
      
      g5 <- invisible(plot_structure(res))
      expect_true(inherits(g5, c("ggplot", "qgraph")))
    }
    
    # --- posterior parameter HDI ---
    g6 <- invisible(plot_parameterHDI(res))
    expect_s3_class(g6, "ggplot")
    
    # --- centrality ---
    g7 <- invisible(plot_centrality(res))
    expect_s3_class(g7, "ggplot")
  }
})




