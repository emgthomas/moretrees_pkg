context("Running moretrees for different input parameter values")
require(moretrees)
#require(testthat)

# Parameter choices to be tested -----------------------------------------------------
params_list <- list(K = 1:2, 
                    m = 0:2,
                    nrestarts = c(1, 2))
params <- do.call(expand.grid, params_list)

# Load data --------------------------------------------------------------------------
data("moretreesExampleEdgelist")
tr <- igraph::graph_from_edgelist(moretreesExampleEdgelist, directed = TRUE)
data("moretreesExampleData")

# Run tests for different parameters -------------------------------------------------

for (i in 1:nrow(params)) {
  K <- params$K[i]
  m <- params$m[i]
  nrestarts <- params$nrestarts[i]
  mod_end <- test_moretrees(moretreesExampleData,
                            tr, K, m, nrestarts)
  test_that(paste0("model output is a list when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_is(mod_end, "moretrees_result")
                     expect_is(mod_end, "list")
                   })
  test_that(paste0("dimension of beta_est is K*3+1 when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_equal(ncol(mod_end$beta_est), K * 3 + 1)
                   })
  test_that(paste0("dimension of beta_moretrees K*3+3 when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_equal(ncol(mod_end$beta_moretrees), K * 3 + 3)
                   })
  test_that(paste0("dimension of beta_ml is K*3+2 when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_equal(ncol(mod_end$beta_ml), K * 3 + 2)
                   })
  test_that(paste0("dimension of theta_est is m*3 when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     if (m == 0) expect_equal(mod_end$theta_est, NULL)
                     if (m > 0) expect_equal(ncol(mod_end$theta_est), m * 3)
                   })
  test_that(paste0("Dimension of theta_ml when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     if (m == 0) expect_equal(mod_end$theta_ml, NULL)
                     if (m > 0) expect_equal(ncol(mod_end$theta_ml), m * 3 + 2)
                   })
  # Test plotting method
  p <- plot(mod_end)
  test_that(paste0("Plotting works when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_is(p, "ggtree")
                     expect_is(p, "ggplot")
                   })
  
  # Test summary methods
  mod_long <- summary(mod_end, compact = FALSE)
  mod_short <- summary(mod_end, compact = TRUE)
  test_that(paste0("Summary works when K = ", K,
                   ", m = ", m, 
                   ", nrestarts = ", nrestarts), {
                     expect_is(mod_long, "summary.moretrees_long")
                     expect_is(mod_short, "summary.moretrees_compact")
                   })
  
  # Test printing methods
  print(mod_long)
  print(mod_short)
  
  # Test printing methods for CLR on MOReTreeS groups
  print(mod_end, coeff_type = "clr", compact = T)
}
