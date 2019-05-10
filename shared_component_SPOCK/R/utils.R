run_sim <- function(n_rep, map_state, n_cov, confounding, alpha, betas, tau, delta, prior, cores = 10){
  ##-- Map ----
  nb_list <- spdep::poly2nb(map_state)
  W <- spdep::nb2mat(nb_list, style = "B")
  n <- length(map_state)

  ##-- Covariates ----
  X <- matrix(rnorm(n = n_cov*n, mean = 0, sd = 1), ncol = n_cov)

  if(confounding){
    X[, n_cov] <- scale(coordinates(map_state)[, 2])
  }

  colnames(X) <- paste0("X", 1:ncol(X))

  ##-- spock ----
  spock_nb <- spock(x = cbind(1, X), map = map_state, use_map = TRUE)
  W_spock <- spdep::nb2mat(spock_nb, style = "B")

  cl <- makeCluster(cores)

  clusterEvalQ(cl = cl, library("INLA"))
  clusterEvalQ(cl = cl, library("spdep"))
  clusterEvalQ(cl = cl, library("sp"))
  clusterEvalQ(cl = cl, library("mvtnorm"))
  clusterExport(cl = cl, c("ricar", "get_trans_est", "run_shared", "try_run_shared"))

  results <- clusterApplyLB(cl = cl, x = 1:n_rep, fun = try_run_shared,
                            map_state = map_state,
                            X = X, confounding = confounding,
                            W = W, W_spock = W_spock,
                            alpha = alpha, betas = betas, tau = tau, delta = delta,
                            prior = prior)

  stopCluster(cl)

  return(results)
}

try_run_shared <- function(i, map_state, X, confounding, W, W_spock, alpha, betas, tau, delta, prior){
  results <- try(expr = run_shared(i, map_state, X, confounding, W, W_spock, alpha, betas, tau, delta, prior))

  if(class(results) == "try-error"){
    results <- try_run_shared(i, map_state, X, confounding, W, W_spock, alpha, betas, tau, delta, prior)
  }

  return(results)
}

run_shared <- function(i, map_state, X, confounding, W, W_spock, alpha, betas, tau, delta, prior){
  set.seed(i)

  X_aux <- cbind(1, X)

  ##-- Spatial effects ----
  sc <- ricar(W = W, num = colSums(W), sig = 1/tau[1])

  Ps <- sc%*%solve(t(sc)%*%sc)%*%t(sc)
  Ps_ort <- diag(nrow(Ps)) - Ps
  s1 <- Ps_ort%*%ricar(W = W, num = colSums(W), sig = 1/tau[2])

  sc_s1 <- cbind(sc, s1)
  Ps2 <- sc_s1%*%solve(t(sc_s1)%*%sc_s1)%*%t(sc_s1)
  Ps2_ort <- diag(nrow(Ps2)) - Ps2

  s2 <- Ps2_ort%*%ricar(W = W, num = colSums(W), sig = 1/tau[3])

  spatial_1 <- sc*delta + s1
  spatial_2 <- sc/delta + s2

  ##-- beta star SPOCK ----
  ll <- solve(t(X_aux)%*%X_aux)%*%t(X_aux)%*%cbind(1, coordinates(map_state))

  bb_1 <- coefficients(lm(spatial_1 ~ coordinates(map_state)))
  bb_2 <- coefficients(lm(spatial_2 ~ coordinates(map_state)))

  beta_star_1 <- as.numeric(c(alpha[1], betas[1, ]) + ll%*%bb_1)
  beta_star_2 <- as.numeric(c(alpha[2], betas[2, ]) + ll%*%bb_2)

  ##-- SMR, expected value and counts----
  n <- nrow(X)
  n_covs <- length(betas)

  srm_1 <- exp(alpha[1] + X%*%betas[1, ] + spatial_1)
  E_1 <- runif(n = n, min = 10, max = 20)
  Y_1 <- rpois(n = n, lambda = E_1*srm_1)

  srm_2 <- exp(alpha[2] + X%*%betas[2, ] + spatial_2)
  E_2 <- runif(n = n, min = 10, max = 20)
  Y_2 <- rpois(n = n, lambda = E_2*srm_2)

  ##-- Dataset ----
  df_pois <- data.frame(ID = 1:n, Y_1 = Y_1, Y_2 = Y_2, E_1 = E_1, E_2 = E_2, X, s1 = s1, s2 = s2, sc = sc)

  ##-- Model ----
  ##-- + Setup INLA ----
  Y <- df_pois[, c("Y_1", "Y_2")]
  E <- df_pois[, c("E_1", "E_2")]

  inla_list <- list()
  inla_list$Y <- cbind(c(Y[, 1], rep(NA, n)),
                       c(rep(NA, n), Y[, 2]))
  E <- c(E[, 1], E[, 2])

  inla_list$alpha1 <- rep(c(1, NA), each = n)
  inla_list$alpha2 <- rep(c(NA, 1), each = n)

  inla_list$psi_gamma <- c(df_pois$ID, rep(NA, n))
  inla_list$psi <- c(rep(NA, n), df_pois$ID)

  inla_list$phi1 <- inla_list$psi_gamma
  inla_list$phi2 <- inla_list$psi

  inla_list$X <- matrix(NA, nrow = 2*n, ncol = n_covs)
  inla_list$X[1:n, 1:(n_covs/2)] <- as.matrix(X)
  inla_list$X[(n+1):(2*n), ((n_covs/2)+1):n_covs] <- as.matrix(X)
  colnames(inla_list$X) <- c(paste0(colnames(X), 1), paste0(colnames(X), 2))

  ##-- Params
  params <- data.frame(param = c("alpha1", "alpha2",
                                 paste0("X", 1:(n_covs/2), "1"), paste0("X", 1:(n_covs/2), "2"),
                                 "Delta", "Beta for psi_gamma", "Precision for psi", "Precision for psi adj", "Precision for phi1", "Precision for phi2"),
                       real = c(alpha,
                                betas[1, ], betas[2, ],
                                delta, delta^2, tau[1]*delta^2, tau[1], tau[2], tau[3]),
                       stringsAsFactors = FALSE)

  ##-- + Fit shared ----
  f_s <- Y ~ -1 + alpha1 + alpha2 + X +
    f(psi, model = "besag", graph = W, hyper = list(theta = list(initial = 0))) +
    f(psi_gamma, copy = "psi", fixed = FALSE, range = c(0, Inf)) +
    f(phi1, model = "besag", graph = W, hyper = list(theta = list(initial = 0, param = prior))) +
    f(phi2, model = "besag", graph = W, hyper = list(theta = list(initial = 0, param = prior)))

  mod <- inla(formula = f_s,
              family = c("poisson", "poisson"),
              data = inla_list,
              E = as.vector(E),
              control.compute = list(waic = T))

  fixed <- mod$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  fixed$beta_star <- c(alpha[1], alpha[2], betas[1,], betas[2,])
  hyperpar <- mod$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.975quant")]
  waic <- data.frame(waic = mod$waic$waic)
  rownames(waic) <- "waic"

  delta_df <- get_trans_est(fun = sqrt, marg = mod$marginals.hyperpar$`Beta for psi_gamma`, name = "Delta", trunc = TRUE)
  tau_psi_df <- get_trans_est(fun = function(x) x/hyperpar$mean[4], marg = mod$marginals.hyperpar$`Precision for psi`, name = "Precision for psi adj", trunc = FALSE)
  hyperpar <- rbind.data.frame(hyperpar, delta_df, tau_psi_df)

  result_list <- lapply(list(fixed, hyperpar, waic),
                        function(x) cbind(seed = seed, confounding = confounding, model = "shared", param = rownames(x), delta = delta, x, stringsAsFactors = FALSE))
  names(result_list) <- c("fixed", "hyperpar", "waic")
  result_list$fixed <- merge(result_list$fixed, params, by = "param")
  result_list$hyperpar <- merge(result_list$hyperpar, params, by = "param")

  ##-- + Fit spock ----
  f_spock <- Y ~ -1 + alpha1 + alpha2 + X +
    f(psi, model = "besag", graph = W_spock, hyper = list(theta = list(initial = 0))) +
    f(psi_gamma, copy = "psi", fixed = FALSE, range = c(0, Inf)) +
    f(phi1, model = "besag", graph = W_spock, hyper = list(theta = list(initial = 0, param = prior))) +
    f(phi2, model = "besag", graph = W_spock, hyper = list(theta = list(initial = 0, param = prior)))

  mod_spock <- inla(formula = f_spock,
                    family = c("poisson", "poisson"),
                    data = inla_list,
                    E = as.vector(E),
                    control.compute = list(waic = T))

  fixed_s <- mod_spock$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  fixed_s$beta_star <- c(beta_star_1[1], beta_star_2[1], beta_star_1[-1], beta_star_2[-1])
  hyperpar_s <- mod_spock$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.975quant")]
  waic_s <- data.frame(waic = mod_spock$waic$waic)
  rownames(waic_s) <- "waic"

  delta_df_s <- get_trans_est(fun = sqrt, marg = mod_spock$marginals.hyperpar$`Beta for psi_gamma`, name = "Delta", trunc = TRUE)
  tau_psi_df_s <- get_trans_est(fun = function(x) x/hyperpar_s$mean[4], marg = mod_spock$marginals.hyperpar$`Precision for psi`, name = "Precision for psi adj", trunc = FALSE)
  hyperpar_s <- rbind.data.frame(hyperpar_s, delta_df_s, tau_psi_df_s)

  result_list_s <- lapply(list(fixed_s, hyperpar_s, waic_s),
                          function(x) cbind.data.frame(seed = seed, confounding = confounding, model = "spock", param = rownames(x), delta = delta, x, stringsAsFactors = FALSE))
  names(result_list_s) <- c("fixed", "hyperpar", "waic")
  result_list_s$fixed <- merge(result_list_s$fixed, params, by = "param")
  result_list_s$hyperpar <- merge(result_list_s$hyperpar, params, by = "param")

  ##-- + Results ----
  results <- Map(rbind, result_list_s, result_list)
  return(results)
}

spock <- function(x, map, coords, map_nb, use_map = !is.null(map)){
  X <- x
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  n <- ncol(H)

  if(use_map) eigenvectors <- (diag(n) - H)%*%coordinates(map)
  else eigenvectors <- (diag(n) - H)%*%coords

  if(use_map) map_nb <- spdep::poly2nb(map)
  n_neigh <- sapply(map_nb, length)

  dist_mat <- as.matrix(dist(eigenvectors))

  nearest_nb <- apply(dist_mat, 2, function(x) order(x, decreasing = FALSE))
  k_nearest_nb <- sapply(1:length(map_nb), function(x) c(as.integer(nearest_nb[1:(n_neigh[x]+1), x]))[-1])
  class(k_nearest_nb) <- "nb"

  A <- spdep::nb2mat(k_nearest_nb, style="B", zero.policy = TRUE)
  A <- A + t(A)
  A[which(A>0, arr.ind=T)] <- 1

  mapa_nb <- apply(A, 1, function(x) which(x > 0))
  class(mapa_nb) <- "nb"

  return(mapa_nb)
}

ricar <- function(W, num, sig = 1, P = NULL){
  n <- ncol(W)
  Q <- -W
  diag(Q) <- num

  if(!is.null(P)) Q <- P%*%Q%*%t(P)

  Q_aux <- eigen(Q)$vectors[, order(eigen(Q)$values)]

  D_aux <- sort(eigen(Q)$values)

  rnd <- rnorm(n-1, 0, sqrt(sig*(1/D_aux[-1])))
  rnd <- Q_aux%*%c(0, rnd)

  return(as.vector(rnd))
}

gamma_prior <- function(mean, var, a = NULL, b = NULL){
  if(all(!is.null(c(a, b)))){
    mean <- a/b
    var <- a/b^2

    return(c(mean = mean, var = var))
  } else{
    a <- (mean^2)/var
    b <- mean/var

    return(c(a = a, b = b))
  }
}

get_trans_est <- function(marg, fun, name, trunc = FALSE, method = "linear"){

  if(trunc){
    marg <- marg[marg[, 1] > 0,]
  }

  delta_marg <- inla.tmarginal(fun = fun, marginal = marg, method = method)
  delta_stat <- inla.zmarginal(delta_marg, silent = T)

  delta_df <- data.frame(mean = delta_stat$mean, sd = delta_stat$sd,
                         `0.025quant` = delta_stat$`quant0.025`, `0.975quant` = delta_stat$`quant0.975`,
                         check.names = FALSE, stringsAsFactors = FALSE)
  rownames(delta_df) <- name
  return(delta_df)
}

extract_stat <- function(m, name){
  fixed <- m$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  fixed <- data.frame(param = rownames(fixed), fixed, stringsAsFactors = FALSE, check.names = FALSE)
  hyperpar <- m$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.975quant")]
  hyperpar <- data.frame(param = rownames(hyperpar), hyperpar, stringsAsFactors = FALSE, check.names = FALSE)
  waic <- data.frame(param = "WAIC", mean = m$waic$waic, sd = NA, `0.025quant` = NA, `0.975quant` = NA, stringsAsFactors = FALSE, check.names = FALSE)
  
  if(any(grepl(x = hyperpar$param, pattern = "Beta"))){
    marg <- get_trans_est(marg = m$marginals.hyperpar$`Beta for phi_beta`, fun = sqrt, name = "Delta", trunc = FALSE, method = "linear")
    hyperpar <- bind_rows(hyperpar, cbind.data.frame(param = "Delta", marg, stringsAsFactors = FALSE))
  }
  
  df_out <- data.frame(model = name, bind_rows(fixed, hyperpar, waic), stringsAsFactors = FALSE, check.names = FALSE) %>%
    rename(lim_inf = `0.025quant`, lim_sup = `0.975quant`) %>%
    mutate(HPD = paste0("[", formatC(x = lim_inf, digits = 2, format = "f"), "; ", formatC(x = lim_sup, digits = 2, format = "f"), "]")) %>%
    select(-lim_inf, -lim_sup) %>%
    mutate(HPD = if_else(param == "WAIC", NA_character_, HPD))
  colnames(df_out) <- c("model", "param", "mean", "sd", "HPD")
  
  return(df_out)
}