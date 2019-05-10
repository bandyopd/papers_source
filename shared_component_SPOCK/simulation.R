##-- Packages ----
library(dplyr)
library(sp)
library(sf)
library(spdep)
library(tigris)
library(mvtnorm)
library(INLA)
library(parallel)

##-- Source ----
source('R/utils.R')

##-- Definitions ----
set.seed(123456)
state_fips <- "06"

n_cov <- 2
alpha <- c(0.5, 0.1)
betas_1 <- c(-0.5, -0.2)
betas_2 <- c(-0.8, -0.4)
betas <- matrix(c(betas_1, betas_2), nrow = 2, byrow = TRUE)
tau <- c(1, 5000, 5000)
delta <- c(1, 1.5, 1.75)

n_rep <- 1000
prior <- gamma_prior(mean = 1000, var = 2e+06)

map_state <- tigris::counties(state = "CA")
map_state@data$region_id <- row.names(map_state@data)
map_state@data$region_inla <- 1:nrow(map_state@data)

results <- list()
results_conf <- list()

for(i in 1:length(delta)){
  t1 <- system.time(
    results[[i]] <- run_sim(n_rep = n_rep, cores = 10,
                            map_state = map_state, n_cov = n_cov, confounding = FALSE,
                            alpha = alpha, betas = betas, tau = tau, delta = delta[i],
                            prior = prior)
  )

  t2 <- system.time(
    results_conf[[i]] <- run_sim(n_rep = n_rep, cores = 10,
                                 map_state = map_state, n_cov = n_cov, confounding = TRUE,
                                 alpha = alpha, betas = betas, tau = tau, delta = delta[i],
                                 prior = prior)
  )

  gc(reset = TRUE)
}

#save(results, results_conf, file = "outputs/results.RData")
#load(file = "outputs/results.RData")

##-- Gathering results ----
results_fixed <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results[[i]], FUN = "[", "fixed") %>% bind_rows())
results_fixed_c <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results_conf[[i]], FUN = "[", "fixed") %>% bind_rows())
results_fixed_j <- bind_rows(results_fixed, results_fixed_c)

results_hyperpar <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results[[i]], FUN = "[", "hyperpar") %>% bind_rows())
results_hyperpar_c <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results_conf[[i]], FUN = "[", "hyperpar") %>% bind_rows())
results_hyperpar_j <- bind_rows(results_hyperpar, results_hyperpar_c)

results_waic <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results[[i]], FUN = "[", "waic") %>% bind_rows())
results_waic_c <- lapply(X = 1:length(delta), FUN = function(i) sapply(X = results_conf[[i]], FUN = "[", "waic") %>% bind_rows())
results_waic_j <- bind_rows(results_waic, results_waic_c)

##-- Coverage ----
fixed_cov <- results_fixed_j %>%
  mutate(cov = if_else(beta_star > `0.025quant` & beta_star < `0.975quant`, 1, 0)) %>%
  group_by(confounding, model, param, delta) %>%
  summarise(cov = mean(cov))

hyperpar_cov <- results_hyperpar_j %>%
  mutate(cov = if_else(real > `0.025quant` & real < `0.975quant`, 1, 0)) %>%
  group_by(confounding, model, param, delta) %>%
  summarise(cov = mean(cov))

waic <- results_waic_j %>%
  reshape2::dcast(seed + confounding + delta ~ model, value.var = "waic", fun.aggregate = sum)
  tidyr::spread(key = model, value = waic) %>%
  mutate(win = if_else(shared > spock, "spock", "shared")) %>%
  group_by(confounding, delta) %>%
  summarise(win_spock = sum(win == "spock"),
            win_shared = sum(win == "shared"))