
### how do i vary the versions of bgms with easybgm 
##### CROSS-SECTIONAL
###-------------
### Estimation checks 
###-------------

test_that("easybgm returns expected structure across valid type–package combos", {
  set.seed(123)
  
  # Subsample small data to stay fast on CRAN
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p <- ncol(dat)
  itr <- 10
  
  if(packageVersion("bgms") >= "0.2.0.0"){
  # Test only core combinations
  combos <- list(
    ### BGGM
    list(type = "continuous", pkg = "BGGM", sv = F, cnt = F),
    list(type = "continuous", pkg = "BGGM", sv = T, cnt = T),
    list(type = "mixed", pkg = "BGGM", sv = T, cnt = T),
    ### BDGRAPH
    list(type = "mixed",      pkg = "BDgraph", sv = F, cnt = F),
    list(type = "continuous",  pkg = "BDgraph", sv = F, cnt = F),
    ### bgms
    list(type = "binary",     pkg = "bgms", sv = F, cnt = F),
    list(type = "binary",     pkg = "bgms", sv = T, cnt = T),
    list(type = "binary",     pkg = "bgms", sv = F, cnt = T),
    list(type = "blume-capel", pkg = "bgms", sv = T, cnt = T),
    list(type = "binary", pkg = "bgms", sv = T, cnt = T, sbm = "Stochastic-Block"),
    list(type = "continuous", pkg = "bgms", sv = F, cnt = F),
    list(type = "continuous", pkg = "bgms", sv = T, cnt = T),
    list(type = "mixed",      pkg = "bgms", sv = T, cnt = T),
    list(type = c("ordinal", "ordinal", "continuous", "continuous", "ordinal"), pkg = "bgms", sv = F, cnt = F)
  )
  } else {
    # Test only core combinations
    combos <- list(
      ### BGGM
      list(type = "continuous", pkg = "BGGM", sv = F, cnt = F),
      list(type = "continuous", pkg = "BGGM", sv = T, cnt = T),
      list(type = "mixed", pkg = "BGGM", sv = T, cnt = T),
      ### BDGRAPH
      list(type = "mixed",      pkg = "BDgraph", sv = F, cnt = F),
      list(type = "continuous",  pkg = "BDgraph", sv = F, cnt = F),
      ### bgms
      list(type = "binary",     pkg = "bgms", sv = F, cnt = F),
      list(type = "binary",     pkg = "bgms", sv = T, cnt = T),
      list(type = "binary",     pkg = "bgms", sv = F, cnt = T),
      list(type = "blume-capel", pkg = "bgms", sv = T, cnt = T),
      list(type = "binary", pkg = "bgms", sv = T, cnt = T, sbm = "Stochastic-Block")
    )
  }
  
  for (cmb in combos) {
    t <- cmb$type
    pkg <- cmb$pkg
    sv <- cmb$sv
    cnt <- cmb$cnt
    if(!is.null(cmb$sbm)) {sbm <- cmb$sbm}
    
    not_cont <- if (length(t) == 1 && t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL

    # bgms defaults to warmup = 2000, which dominates the runtime at these tiny
    # iteration counts. 300 is the shortest warmup bgms does not warn about.
    # BGGM and BDgraph have no warmup argument.
    # bgms also defaults to cores = parallel::detectCores(); CRAN caps check
    # processes at 2 cores, so pin it rather than claiming every core on the
    # machine. detectCores() can also return NA on restricted hosts, which
    # would otherwise make the fit fail.
    extra <- if (identical(pkg, "bgms")) list(warmup = 300, cores = 2L) else list()

    base_args <- list(
      data       = dat,
      type       = t,
      package    = pkg,
      iter       = itr,          # tiny for speed
      save       = sv,
      centrality = cnt,
      progress   = FALSE
    )

    if(length(t) == 1 && t == "blume-capel"){
      suppressWarnings({
        res <- do.call(easybgm, c(base_args, extra,
                                  list(not_cont = not_cont,
                                       baseline_category = 2)))
      })} else if(!is.null(cmb$sbm)){
        suppressWarnings({
          res <- do.call(easybgm, c(base_args, extra,
                                    list(edge_prior = sbm)))
        })
      } else {
        suppressWarnings({
          res <- do.call(easybgm, c(base_args, extra,
                                    list(not_cont = not_cont)))
        })
      }
    
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
    
    expect_no_error(summary(res))
    
    if(sv == TRUE && pkg == "BGGM") {
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(itr, k))
      expect_equal(dim(res$centrality),  c(itr, p))
    } 
    if(cnt == TRUE && pkg == "bgms"){
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(4*itr, k))
      expect_equal(dim(res$centrality),  c(4*itr, p))
    }
    if(!is.null(cmb$sbm)){
      expect_equal(length(res$sbm), 4)
    }
    print(paste0("Finished easybgm: Package: ", cmb$pkg, "; Type: ", cmb$type, "; Centrality: ", cmb$cnt))
    
  }
})

###-------------
### Plotting functions test
###-------------

test_that("plotting functions work across valid type–package combos", {
  set.seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p   <- ncol(dat)
  
  combos <- list(
    list(type = "continuous", pkg = "BGGM"),
   # list(type = "mixed",      pkg = "BDgraph"),
    list(type = "binary",     pkg = "bgms")
  )
  
  for (cmb in combos) {
    t   <- cmb$type
    pkg <- cmb$pkg
    not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL
    
    
    if(pkg == "BDgraph") {
      suppressMessages({
        res <- easybgm(
          data       = dat,
          type       = t,
          package    = pkg,
          iter       = 10,
          save       = FALSE,
          centrality = TRUE,
          progress   = FALSE,
          not_cont   = not_cont
        )
      }) 
    } else {
      # bgms defaults to warmup = 2000 and cores = detectCores(); BGGM has
      # neither argument. Cores are pinned to 2 for CRAN's check limit.
      extra <- if (identical(pkg, "bgms")) list(warmup = 300, cores = 2L) else list()
      suppressWarnings({
        res <- do.call(easybgm, c(
          list(
            data       = dat,
            type       = t,
            package    = pkg,
            iter       = 100,
            save       = TRUE,
            centrality = TRUE,
            progress   = FALSE,
            not_cont   = not_cont
          ), extra))
      })
    }
    
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
    if(pkg != "BDgraph"){
      g6 <-    suppressWarnings({invisible(plot_parameterHDI(res))})
      expect_s3_class(g6, "ggplot")
      
      # --- centrality ---
      g7 <- invisible(plot_centrality(res))
      expect_s3_class(g7, "ggplot")
    }
  }
})

# # TEst only possible to include post 0.2.0.0 version
# test_that("easybgm defaults to bgms for all data types", {
#   data("Wenchuan", package = "bgms")
#   dat <- na.omit(Wenchuan)[1:20, 1:5]
#   
#   suppressWarnings({
#     res <- easybgm(dat, type = "continuous", iter = 10, progress = FALSE)
#   })
#   expect_true("package_bgms" %in% class(res))
#   
#   suppressWarnings({
#     res2 <- easybgm(dat, type = "ordinal", iter = 10, progress = FALSE)
#   })
#   expect_true("package_bgms" %in% class(res2))
# })

test_that("easybgm errors for BDgraph continuous with missing data", {
  data("Wenchuan", package = "bgms")
  dat_with_na <- Wenchuan[1:20, 1:5]  # Wenchuan has NAs
  
  expect_error(
    easybgm(dat_with_na, type = "continuous", package = "BDgraph",
            iter = 10, progress = FALSE),
    "missing values"
  )
})


##### NETWORK COMPARISON

# test_that("easybgm_compare errors for continuous/mixed without BGGM", {
#   data("Wenchuan", package = "bgms")
#   dat <- na.omit(Wenchuan)[1:20, 1:5]
#   group_dat <- list(dat[1:10, ], dat[11:20, ])
#   expect_error(
#     suppressWarnings(
#       easybgm_compare(group_dat, type = "continuous", package = "bgms")
#     ),
#     "bgms can only fit 'binary'"
#   )
#   expect_error(
#     suppressWarnings(
#       easybgm_compare(group_dat, type = "mixed", package = "bgms")
#     ),
#     "bgms can only fit 'binary'"
#   )
# })


test_that("easybgm_compare accepts a per-variable type vector", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:30, 1:3]
  grp <- rep(1:2, length.out = nrow(dat))
  fit <- suppressWarnings(
    easybgm_compare(dat, type = rep("ordinal", 3), group_indicator = grp,
                    iter = 50, warmup = 300, cores = 2L, progress = FALSE)
  )
  expect_s3_class(fit, "package_bgms_compare")

  # a vector whose length does not match the number of columns is rejected
  expect_error(
    suppressWarnings(
      easybgm_compare(dat, type = rep("ordinal", 2), group_indicator = grp,
                      iter = 50, progress = FALSE)
    ),
    "must equal"
  )
})


test_that("easybgm_compare returns expected structure across valid type–package combos", {
  set.seed(123)
  
  # Subsample small data to stay fast on CRAN
  data("Wenchuan", package = "bgms")
  dat <- as.data.frame(na.omit(Wenchuan)[1:90, 1:5])
  p <- ncol(dat)
  itr <- 10
  
  # Test only core combinations
  combos <- list(
    ### BGGM
    list(type = "continuous", pkg = "BGGM", sv = F),
    list(type = "continuous", pkg = "BGGM", sv = T),
    list(type = "mixed", pkg = "BGGM", sv = T),
    ### bgms
    list(type = "binary",     pkg = "bgms", sv = F),
    list(type = "binary",     pkg = "bgms", sv = T),
    list(type = "binary",     pkg = "bgms", sv = T, multi_group = T)
  )
  
  for (cmb in combos) {
    t <- cmb$type
    pkg <- cmb$pkg
    sv <- cmb$sv
    
    # bgms defaults to warmup = 2000 and cores = detectCores(); BGGM has
    # neither argument. Cores are pinned to 2 for CRAN's check limit.
    extra <- if (identical(pkg, "bgms")) list(warmup = 300, cores = 2L) else list()

    if(!is.null(cmb$multi_group)){
      group <- rep(c(1, 2, 3), each = 30)

      suppressMessages({
        res <- do.call(easybgm_compare, c(
          list(
            data       = dat,
            type       = t,
            package    = pkg,
            iter       = itr,          # tiny for speed
            save       = sv,
            group_indicator = group,
            progress   = FALSE
          ), extra))
      })
    } else {
      group_dat <- list(dat[1:45, ], dat[46:90, ])
      not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL

      suppressWarnings({
        res <- do.call(easybgm_compare, c(
          list(
            data       = group_dat,
            type       = t,
            package    = pkg,
            iter       = itr,          # tiny for speed
            save       = sv,
            progress   = FALSE,
            not_cont   = not_cont
          ), extra))
      })
    }
    # --- class check ---
    expect_true(inherits(res, c("easybgm_compare")))
    expect_true(any(grepl("package_", class(res))))  # backend tag present
    
    # --- field presence check ---
    # multi-group bgms comparisons return pairwise group differences instead
    # of a single difference matrix
    par_field <- if (is.null(cmb$multi_group)) "parameters" else "pairwise_group_differences"
    expect_true(all(c(par_field, "inc_probs", "inc_BF", "structure", "model") %in% names(res)))
    
    # --- dimensions check ---
    if (is.null(cmb$multi_group)) expect_equal(dim(res$parameters), c(p, p))
    expect_equal(dim(res$inc_probs),  c(p, p))
    expect_equal(dim(res$inc_BF),     c(p, p))
    expect_equal(dim(res$structure),  c(p, p))
    
    # --- sanity check ---
    expect_false(all(is.na(res[[par_field]])))
    expect_false(all(is.na(res$inc_probs)))
    
    if(sv == TRUE && pkg != "bgms") {
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(itr, k))
    }
    if(sv == TRUE && pkg == "bgms"){
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(4*itr, k))
    }
    
    
    print(paste0("Finished easybgm_compare: Package: ", cmb$pkg, "; Type: ", cmb$type))
  }
})


###-------------
### Regression checks for bgms >= 0.2.0.0 extraction
###-------------

test_that("bgms centrality uses the lower-triangle edge ordering", {
  # bgms has stored pairwise interactions in lower-triangle order since at least
  # 0.1.6.3, so this holds on both supported bgms versions.
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]

  res <- suppressWarnings(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300, cores = 2L,
            save = TRUE, centrality = TRUE, progress = FALSE)
  )
  p <- ncol(res$parameters)

  # bgms stores pairwise interactions in lower-triangle column order, the same
  # order res$parameters is filled from.
  expected <- t(apply(res$samples_posterior, 1, function(r)
    rowSums(abs(vector2matrix(r, p, bycolumn = FALSE)))))
  expect_equal(unname(res$centrality), unname(expected))

  # The upper-triangle fill (BGGM's order) genuinely differs, so the check above
  # would fail if the ordering regressed.
  wrong <- t(apply(res$samples_posterior, 1, function(r)
    rowSums(abs(vector2matrix(r, p, bycolumn = TRUE)))))
  expect_false(isTRUE(all.equal(unname(expected), unname(wrong))))
})

test_that("structure is the median probability model when save = FALSE", {
  # both supported bgms versions report inclusion probabilities the same way
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:4]

  res <- suppressWarnings(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300, cores = 2L,
            save = FALSE, progress = FALSE)
  )
  expect_equal(unname(res$structure), unname(1 * (res$inc_probs > 0.5)))
  expect_true(all(diag(res$structure) == 0))
})

test_that("entry points accept raw bgms fit objects", {
  # Exercised on both supported bgms versions: an S7 object on bgms >= 0.2.0.0
  # and a plain S3 list on 0.1.6.3. The prior constructors only exist on the
  # newer version, so the SBM prior is specified in whichever form applies.
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:4]

  fit <- suppressWarnings(
    bgms::bgm(dat, iter = 50, warmup = 300, chains = 2, cores = 2L,
              display_progress = FALSE)
  )
  # these two list methods used to fail on the missing `save` fit argument
  expect_no_error(suppressWarnings(plot_centrality(list(fit, fit))))
  expect_no_error(suppressWarnings(plot_prior_sensitivity(list(fit, fit))))

  # clusterBayesfactor used to call names()<- on the fit, which an S7 object
  # does not allow, and read $sbm, which a raw fit does not carry
  sbm_arg <- if (packageVersion("bgms") >= "0.2.0.0") {
    bgms::sbm_prior()
  } else {
    "Stochastic-Block"
  }
  fit_sbm <- suppressWarnings(
    bgms::bgm(dat, edge_prior = sbm_arg, iter = 50, warmup = 300,
              chains = 2, cores = 2L, display_progress = FALSE)
  )
  expect_no_error(suppressWarnings(clusterBayesfactor(fit_sbm)))
})

test_that("legacy interaction_scale does not raise a bgms deprecation warning", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:3]

  w <- character(0)
  withCallingHandlers(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300, cores = 2L,
            progress = FALSE, interaction_scale = 2.5),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("deprecat", w, ignore.case = TRUE)))
})

###-------------
### Blume-Capel main effects
###-------------

test_that("easybgm reports Blume-Capel linear and quadratic effects", {
  skip_if_not(packageVersion("bgms") >= "0.2.0.0")

  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:60, 1:4]
  p <- ncol(dat)
  itr <- 20

  res <- suppressWarnings(easybgm(dat, type = "blume-capel", package = "bgms",
                                  baseline_category = 3, iter = itr,
                                  warmup = 300, cores = 2L,
                                  progress = FALSE, save = TRUE))

  bc <- res$blume_capel_parameters
  expect_s3_class(bc, "data.frame")
  # two parameters, linear and quadratic, for every variable
  expect_equal(nrow(bc), 2 * p)
  expect_equal(unique(bc$Variable), colnames(dat))
  expect_equal(unique(bc$Effect), c("linear", "quadratic"))
  expect_false(any(is.na(bc$Estimate)))

  # The baseline is reported on the scale of the input data. bgms recodes
  # scores to start at 0 and shifts the baseline with them, so reading
  # args$baseline_category without adding back the shift understates it.
  expect_true(all(bc[["Baseline Category"]] == 3))

  # the table summarises exactly the draws that were saved
  expect_equal(unname(colMeans(res$samples_blume_capel)), bc$Estimate)
  expect_equal(colnames(res$samples_blume_capel),
               paste0(bc$Variable, " (", bc$Effect, ")"))
  expect_equal(nrow(res$samples_blume_capel), 4 * itr)

  # credible interval brackets the posterior mean
  expect_true(all(bc[["Lower 2.5%"]] <= bc$Estimate))
  expect_true(all(bc$Estimate <= bc[["Upper 97.5%"]]))

  # $thresholds holds the same two parameters, now named for what they are
  expect_equal(colnames(res$thresholds), c("linear", "quadratic"))
  expect_equal(as.vector(t(res$thresholds[unique(bc$Variable), ])), bc$Estimate)

  # summary carries the table and print emits it exactly once
  s <- summary(res)
  expect_s3_class(s$blume_capel_parameters, "data.frame")
  out <- capture.output(print(s))
  expect_equal(sum(grepl("BLUME-CAPEL MAIN EFFECTS", out)), 1)
  expect_equal(sum(grepl("BLUME-CAPEL MAIN EFFECTS", capture.output(print(res)))), 1)
})

test_that("Blume-Capel effects are correct when mixed with other variable types", {
  skip_if_not(packageVersion("bgms") >= "0.2.0.0")

  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:60, 1:4]
  # a continuous variable makes the bgms argument vectors shorter than p:
  # baseline_category and blume_capel_shift are indexed over discrete
  # variables only, so a positional read misaligns them.
  type <- c("continuous", "blume-capel", "ordinal", "blume-capel")

  res <- suppressWarnings(easybgm(dat, type = type, package = "bgms",
                                  baseline_category = 4, iter = 20,
                                  warmup = 300, cores = 2L, progress = FALSE))

  bc <- res$blume_capel_parameters
  expect_equal(unique(bc$Variable), colnames(dat)[type == "blume-capel"])
  expect_equal(nrow(bc), 4)
  expect_true(all(bc[["Baseline Category"]] == 4))

  # Blume-Capel and ordinal rows share one matrix, so the column headers
  # cannot describe both and the per-row meaning is recorded instead
  disc <- res$thresholds$discrete
  expect_equal(attr(disc, "variable_type"),
               c("blume-capel", "ordinal", "blume-capel"))
  expect_equal(as.vector(t(disc[unique(bc$Variable), 1:2])), bc$Estimate)
})

test_that("no Blume-Capel output when no variable is Blume-Capel", {
  skip_if_not(packageVersion("bgms") >= "0.2.0.0")

  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:60, 1:4]

  res <- suppressWarnings(easybgm(dat, type = "ordinal", package = "bgms",
                                  iter = 20, warmup = 300, cores = 2L,
                                  progress = FALSE, save = TRUE))

  expect_null(res$blume_capel_parameters)
  expect_null(res$samples_blume_capel)
  expect_null(attr(res$thresholds, "variable_type"))
  expect_null(summary(res)$blume_capel_parameters)
  expect_equal(sum(grepl("BLUME-CAPEL", capture.output(print(res)))), 0)
})

test_that("Blume-Capel samples are only stored when save = TRUE", {
  skip_if_not(packageVersion("bgms") >= "0.2.0.0")

  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:60, 1:4]

  res <- suppressWarnings(easybgm(dat, type = "blume-capel", package = "bgms",
                                  baseline_category = 2, iter = 20,
                                  warmup = 300, cores = 2L,
                                  progress = FALSE, save = FALSE))

  expect_null(res$samples_blume_capel)
  expect_s3_class(res$blume_capel_parameters, "data.frame")
  # the interval does not depend on save = TRUE
  expect_false(any(is.na(res$blume_capel_parameters[["Lower 2.5%"]])))
})


###-------------
### Regression checks against the raw bgms fit (easybgm 0.5.1)
###-------------

test_that("mixed-type pairwise estimates are placed on the edge they belong to", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  set.seed(1)
  n <- 100
  dat <- data.frame(A = sample(0:2, n, TRUE), B = rnorm(n),
                    C = sample(0:2, n, TRUE), D = rnorm(n))

  res <- suppressWarnings(
    easybgm(dat, type = c("ordinal", "continuous", "ordinal", "continuous"),
            iter = 50, warmup = 300, cores = 2L, save = TRUE, centrality = TRUE,
            progress = FALSE)
  )
  # the returned object carries the bgms fit it was built from
  expect_true(inherits(res$packagefit, "bgms"))
  draws <- bgms::extract_pairwise_interactions(res$packagefit)
  means <- colMeans(draws)

  # bgms groups the pairs by variable type, so a positional fill would misplace them
  vars <- colnames(dat)
  expect_false(identical(names(means), combn(vars, 2, paste, collapse = "-")))

  for (i in 1:3) for (j in (i + 1):4) {
    pair <- intersect(c(paste(vars[i], vars[j], sep = "-"),
                        paste(vars[j], vars[i], sep = "-")), names(means))
    expect_equal(res$parameters[i, j], means[[pair]])
    expect_equal(res$parameters[j, i], means[[pair]])
  }

  # strength centrality places each draw on its edge by name
  expected <- t(apply(draws, 1, function(r) rowSums(abs(vector2matrix_named(r, vars)))))
  expect_equal(unname(res$centrality), unname(expected))

  # summary shows every edge its own R-hat
  rhat <- bgms::extract_rhat(res$packagefit)$pairwise
  s <- summary(res)$parameters
  own <- vapply(strsplit(s$Relation, "-"), function(v)
    rhat[[intersect(c(paste(v, collapse = "-"), paste(rev(v), collapse = "-")), names(rhat))]],
    numeric(1))
  expect_equal(s[[grep("^Convergence", colnames(s), value = TRUE)]], round(own, 3))

  # each Monte Carlo interval of an inclusion BF uses that edge's own MCSE
  ind <- res$packagefit$posterior_summary_indicator
  pairs <- which(lower.tri(res$inc_BF), arr.ind = TRUE)
  forward <- paste(vars[pairs[, "col"]], vars[pairs[, "row"]], sep = "-")
  reversed <- paste(vars[pairs[, "row"]], vars[pairs[, "col"]], sep = "-")
  own_mcse <- ind$mcse[ifelse(forward %in% rownames(ind), match(forward, rownames(ind)),
                              match(reversed, rownames(ind)))]
  bf <- res$inc_BF[lower.tri(res$inc_BF)]
  p <- res$inc_probs[lower.tri(res$inc_probs)]
  se_log <- own_mcse / (p * (1 - p))
  upper <- ifelse(is.finite(se_log), exp(log(bf) + stats::qnorm(0.975) * se_log), NA_real_)
  expect_false(all(is.na(upper)))
  expect_equal(res$MCSE_BF$upper, upper)
  expect_equal(rownames(res$MCSE_BF), forward)
})

test_that("vector2matrix_named matches either orientation and falls back by position", {
  expected <- vector2matrix(c(1, 2, 3), p = 3)
  expect_equal(vector2matrix_named(c("C-B" = 3, "A-B" = 1, "A-C" = 2), c("A", "B", "C")),
               expected)
  expect_warning(m <- vector2matrix_named(c(1, 2, 3), c("A", "B", "C")), "by position")
  expect_equal(m, expected)
})

test_that("inclusion Bayes factors match bgms for every edge prior", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  set.seed(2)
  dat <- as.data.frame(matrix(sample(0:2, 100 * 5, TRUE), 100, 5))
  names(dat) <- LETTERS[1:5]

  cases <- list(
    list(prior = bgms::bernoulli_prior(), save = TRUE),
    list(prior = bgms::beta_bernoulli_prior(alpha = 2, beta = 3), save = TRUE),
    list(prior = bgms::sbm_prior(), save = TRUE),
    # unequal within- and between-block shapes: the prior odds cannot be read
    # off the posterior partition; checked on both extraction branches
    list(prior = bgms::sbm_prior(alpha = 9, beta = 1, alpha_between = 1, beta_between = 9),
         save = TRUE),
    list(prior = bgms::sbm_prior(alpha = 9, beta = 1, alpha_between = 1, beta_between = 9),
         save = FALSE)
  )
  for (cs in cases) {
    res <- suppressWarnings(
      easybgm(dat, type = "ordinal", iter = 50, warmup = 300, cores = 2L,
              save = cs$save, progress = FALSE, edge_prior = cs$prior)
    )
    expected <- bgms::extract_inclusion_bf(res$packagefit)
    diag(expected) <- 0
    expect_equal(unname(res$inc_BF), unname(expected))
    # the prior inclusion probability the prior sensitivity plot uses
    prior_probs <- bgms::extract_prior_inclusion_probabilities(res$packagefit)
    expect_equal(res$edge.prior[1], prior_probs[2, 1])
  }
})

test_that("two-group bgms comparison differences are group 2 minus group 1", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  set.seed(4)
  n <- 60
  dat <- as.data.frame(matrix(sample(0:2, 2 * n * 4, TRUE), 2 * n, 4))
  names(dat) <- LETTERS[1:4]

  res <- suppressWarnings(
    easybgm_compare(list(dat[1:n, ], dat[n + 1:n, ]), type = "ordinal",
                    iter = 50, warmup = 300, cores = 2L, progress = FALSE)
  )
  expect_true(inherits(res$packagefit, "bgmCompare"))
  gp <- bgms::extract_group_params(res$packagefit)$pairwise_effects_groups
  # matrix positions of each named pair
  ij <- matrix(match(do.call(rbind, strsplit(rownames(gp), "-")), names(dat)), ncol = 2)

  expect_equal(res$parameters[ij], unname(gp[, 2] - gp[, 1]))
  expect_equal(res$parameters_g1[ij], unname(gp[, 1]))
  expect_equal(res$parameters_g2[ij], unname(gp[, 2]))

  # summary shows every edge its own R-hat
  s <- summary(res)$parameters
  pd <- summary(res$packagefit)$pairwise_diff
  rhat <- stats::setNames(pd$Rhat, sub(" \\(diff[0-9]+\\)$", "", pd$parameter))
  expect_equal(s$Convergence, unname(round(rhat[s$Relation], 3)))

  # the same holds when the two groups are given through group_indicator, which
  # takes the multi-group path with a single contrast per edge
  res_gi <- suppressWarnings(
    easybgm_compare(dat, type = "ordinal", group_indicator = rep(1:2, each = n),
                    iter = 50, warmup = 300, cores = 2L, progress = FALSE)
  )
  gp_gi <- bgms::extract_group_params(res_gi$packagefit)$pairwise_effects_groups
  ij_gi <- matrix(match(do.call(rbind, strsplit(rownames(gp_gi), "-")), names(dat)), ncol = 2)
  expect_equal(res_gi$parameters[ij_gi], unname(gp_gi[, 2] - gp_gi[, 1]))
  expect_true("Average Difference" %in% colnames(summary(res_gi)$parameters))
  expect_no_error(suppressWarnings(plot_network(res_gi)))
})

test_that("multi-group comparisons report pairwise group differences, not averaged contrasts", {
  skip_if(packageVersion("bgms") < "0.2.0.0")
  set.seed(3)
  n <- 60
  dat <- as.data.frame(matrix(sample(0:2, 3 * n * 4, TRUE), 3 * n, 4))
  names(dat) <- LETTERS[1:4]

  res <- suppressWarnings(
    easybgm_compare(dat, type = "ordinal", group_indicator = rep(1:3, each = n),
                    iter = 50, warmup = 300, cores = 2L, progress = FALSE)
  )
  gp <- bgms::extract_group_params(res$packagefit)$pairwise_effects_groups
  d <- res$pairwise_group_differences

  expect_equal(colnames(d), c("group2 - group1", "group3 - group1", "group3 - group2"))
  expect_equal(rownames(d), rownames(gp))
  expect_equal(unname(d[, "group2 - group1"]), unname(gp[, 2] - gp[, 1]))
  expect_equal(unname(d[, "group3 - group1"]), unname(gp[, 3] - gp[, 1]))
  expect_equal(unname(d[, "group3 - group2"]), unname(gp[, 3] - gp[, 2]))

  # the across-group estimate is the mean of the group estimates
  ij <- matrix(match(do.call(rbind, strsplit(rownames(gp), "-")), names(dat)), ncol = 2)
  expect_equal(res$overall_estimate[ij], unname(rowMeans(gp)))

  # summary shows every edge its own (baseline) R-hat
  s_rhat <- summary(res)$parameters
  pb <- summary(res$packagefit)$pairwise
  rhat <- stats::setNames(pb$Rhat, pb$parameter)
  expect_equal(s_rhat$Convergence, unname(round(rhat[s_rhat$Relation], 3)))

  # the per-contrast coefficients keep their bgms labels, and no field holds a
  # single difference per edge across the three groups
  expect_true(all(grepl(" \\(diff[0-9]+\\)$", res$contrast_coefficients$parameter)))
  expect_null(res$parameters)
  s <- summary(res)
  expect_false("Average Difference" %in% colnames(s$parameters))
  expect_error(suppressWarnings(plot_network(res)), "pairwise_group_differences")

  # the raw fit is plotted as a two-group comparison, which is flagged
  w <- capture_warnings(plot_network(res$packagefit))
  expect_true(any(grepl("more than two groups", w)))
})

test_that("an explicit package = 'bgms' survives a per-variable type vector", {
  set.seed(1); n <- 40
  dat <- data.frame(A = sample(1:3, n, TRUE), B = sample(1:3, n, TRUE),
                    C = sample(1:3, n, TRUE))
  grp <- rep(1:2, each = n / 2)
  expect_no_warning(
    easybgm_compare(dat, type = c("ordinal", "ordinal", "ordinal"),
                    package = "bgms", group_indicator = grp, iter = 20)
  )
})
