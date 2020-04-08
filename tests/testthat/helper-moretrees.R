test_moretrees <- function(moretreesExampleData, tr,
                           K, m, nrestarts) {
  
  # Load data -------------------------------------------------------------------
  Xcase <- moretreesExampleData$Xcase[ , 1:K, drop = F]
  Xcontrol <- moretreesExampleData$Xcontrol[ , 1:K, drop = F]
  if (m > 0) {
    Wcase <- moretreesExampleData$Wcase[ , 1:m, drop = F]
    Wcontrol <- moretreesExampleData$Wcontrol[ , 1:m, drop = F]
  } else {
    Wcase <- Wcontrol <- NULL
  }
  outcomes <- moretreesExampleData$outcomes
  
  # Run algorithm ----------------------------------------------------------------------
  expect_warning(mod_start <- moretrees(Xcase = Xcase, Xcontrol = Xcontrol,
                                        outcomes = outcomes,
                                        tr = tr,
                                        nrestarts = nrestarts,
                                        get_ml = F,
                                        max_iter = 100,
                                        log_restarts = F))
  # strip out unnecessary parts of initial values
  vi_params_init <- mod_start$mod$vi_params
  vi_params_init[c("delta", "Omega", "Omega_inv", "Omega_det")] <- NULL
  hyperparams_init <- mod_start$mod$hyperparams
  hyperparams_init[c("omega", "ELBO", "g_eta")] <- NULL
  # run next model using initial values from previous model
  log_dir <- "restart_logs"
  dir.create(log_dir)
  expect_message(mod_end <- moretrees(Xcase = Xcase, Xcontrol = Xcontrol,
                                      Wcase = Wcase, Wcontrol = Wcontrol,
                                      outcomes = outcomes,
                                      vi_params_init = vi_params_init,
                                      hyperparams_init = hyperparams_init,
                                      tr = tr,
                                      max_iter = 100,
                                      nrestarts = nrestarts,
                                      get_ml = T,
                                      log_restarts = T,
                                      log_dir = log_dir))
  unlink(log_dir, recursive = T)
  
  return(mod_end)
  
}