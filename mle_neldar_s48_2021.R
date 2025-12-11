rm(list=ls())
library(deSolve)
library(ggplot2)
library(purrr)
library(dplyr)
library(bbmle)

mod4 <- function(times, outputV, parms){
  with(as.list(c(outputV, parms)), {
    arriv_rate = eval(D(expression(L / (1 + exp(-(a + k * x)))), "x"), 
                      list(L = L, a = a, k = k, x = times))
    N = S + I + R
    dS = arriv_rate - beta * S * I / N
    dI = beta * S * I / N - gamma * I
    dR = (1-mu) * gamma * I
    dD = mu * gamma * I
    res <- c(dS, dI, dR, dD)
    return(list(res, pop_size = N))
    
  })
}

load("data/arrival_func_params.RData")
load("data/2021.RData")
morta21 <- morta21 %>% filter(julien_date < 150)
time_vec <- seq(min(morta21$julien_date)-1, 150)
SIRfun_nb <- function(log.beta, log.gamma, log.mu, log.size, log.I0){
  parms21 <- c(beta = exp(log.beta), gamma = exp(log.gamma), mu = plogis(log.mu), 
               L = logistic_func[logistic_func$year == 2021, "L"], 
               a = logistic_func[logistic_func$year == 2021, "a"], 
               k = logistic_func[logistic_func$year == 2021, "k"])
  init_pop_size <- as.numeric(parms21["L"] / (1 + exp(-(parms21["a"] + parms21["k"] * min(time_vec)))))
  init_pop <- c(S = init_pop_size * (1 - plogis(log.I0)),
                I = init_pop_size * plogis(log.I0), R = 0, D = 0)
  ss <- try(as.data.frame(rk4(init_pop, time_vec, mod4, parms21)))
  ss_sub <- ss %>% filter(time %in% morta21$julien_date)
  
  if(inherits(ss_sub, "try-error")){
    return(NA)
  } else {
    return(-sum(dnbinom(morta21$cum_death, mu=ss_sub$D, size=exp(log.size), 
                        log = TRUE)))
  }
}

# initial values
init.params <- expand.grid(
  log.beta = log(seq(0.5, 1.5, by = 0.2)),
  log.gamma = log(1/seq(1, 15, by = 2)),
  log.mu = qlogis(seq(0.29, 0.89, by = 0.2)),
  log.size = log(c(0.1, 1, 50)),
  log.I0 = qlogis(c(0.001, 0.01, 0.1, 0.3, 0.5))
)

mfit_ls = list()
mfit_ls <- map(1:nrow(init.params), function(idx) {
  start_nb <- list(log.beta = init.params$log.beta[idx],
                   log.gamma = init.params$log.gamma[idx],
                   log.mu = init.params$log.mu[idx],
                   log.size = init.params$log.size[idx], 
                   log.I0 = init.params$log.I0[idx])
  mfit <- mle2(SIRfun_nb, start_nb, method = "Nelder-Mead", control = list(maxit=1e6))
  print(paste0("progress: ", 100*idx/nrow(init.params), "%"))
  return(mfit)
})
# save(mfit_ls, file = "data/fits_mle_nelder_s48_2021.RData")
load("data/fits_mle_nelder_s48_2021.RData")

successful_indices <- which(!sapply(mfit_ls, is.character))

if(length(successful_indices)==0){
  print("warning: there is no successful search")
}
start_ls <- lapply(mfit_ls, function (mfit) {
  if(!is.character(mfit)){
    mfit@call$start
  }
})

coef_ls <- lapply(mfit_ls, function(mfit){
  if(!is.character(mfit)){
    mfit@coef
  }
})

cofint_ls <- lapply(mfit_ls, function(mfit){
  if(!is.character(mfit)){
    bbmle::confint(mfit, method = "quad")
  }
})

# Create a combined dataframe
combined_df <- data.frame()

for(i in successful_indices) {
  # Extract starting values
  start_vals <- unlist(start_ls[[i]])
  
  # Extract coefficients
  coef_vals <- coef_ls[[i]]
  
  # Extract confidence intervals
  confint_vals <- cofint_ls[[i]]
  
  # Create row with all information
  row_data <- data.frame(
    index = i,
    # Starting values
    start_log.beta = start_vals["log.beta"],
    start_log.gamma = start_vals["log.gamma"], 
    start_log.mu = start_vals["log.mu"],
    start_log.size = start_vals["log.size"],
    start_log.I0 = start_vals["log.I0"],
    # Fitted coefficients
    coef_log.beta = coef_vals["log.beta"],
    coef_log.gamma = coef_vals["log.gamma"],
    coef_log.mu = coef_vals["log.mu"],
    coef_log.size = coef_vals["log.size"],
    coef_log.I0 = coef_vals["log.I0"],
    # Confidence intervals (lower bounds)
    ci_lower_log.beta = confint_vals["log.beta", 1],
    ci_lower_log.gamma = confint_vals["log.gamma", 1],
    ci_lower_log.mu = confint_vals["log.mu", 1],
    ci_lower_log.size = confint_vals["log.size", 1],
    ci_lower_log.I0 = confint_vals["log.I0", 1],
    # Confidence intervals (upper bounds)
    ci_upper_log.beta = confint_vals["log.beta", 2],
    ci_upper_log.gamma = confint_vals["log.gamma", 2],
    ci_upper_log.mu = confint_vals["log.mu", 2],
    ci_upper_log.size = confint_vals["log.size", 2],
    ci_upper_log.I0 = confint_vals["log.I0", 2],
    # Log-likelihood
    logLik = -mfit_ls[[i]]@min
  )
  
  combined_df <- rbind(combined_df, row_data)
}
# Reset row names
rownames(combined_df) <- NULL

# Transform coefficients to original scale
combined_df_t <- combined_df %>%
  mutate(
    # Transform coefficients to original scale
    beta = exp(coef_log.beta),
    gamma = exp(coef_log.gamma),
    mu = plogis(coef_log.mu),
    size = exp(coef_log.size),
    I0 = plogis(coef_log.I0),
    # Transform confidence intervals
    beta_lower = exp(ci_lower_log.beta),
    beta_upper = exp(ci_upper_log.beta),
    gamma_lower = exp(ci_lower_log.gamma),
    gamma_upper = exp(ci_upper_log.gamma),
    mu_lower = plogis(ci_lower_log.mu),
    mu_upper = plogis(ci_upper_log.mu),
    size_lower = exp(ci_lower_log.size),
    size_upper = exp(ci_upper_log.size),
    I0_lower = plogis(ci_lower_log.I0),
    I0_upper = plogis(ci_upper_log.I0),
    # Calculate infectious period
    infectious_period = 1/gamma,
    inf_period_lower = 1/gamma_upper,
    inf_period_upper = 1/gamma_lower
  )

combined_df_clean <- combined_df_t %>% 
  filter(complete.cases(.))

# Reshape data for plotting multiple parameters
plot_data <- combined_df_clean %>%
  select(beta, infectious_period, mu, size, I0) %>%
  tidyr::pivot_longer(cols = everything(), 
                      names_to = "parameter", 
                      values_to = "value")

ggplot(plot_data, aes(x = parameter, y = value)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Maximum likelihood estimation by Nelder-Mead algorithm\nfrom different initial values (from day 46)") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~parameter, scales = "free", ncol = 2) +
  theme(axis.text.x = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_blank())
ggsave("figures/mle_nel_2022_start46.png", width = 5, height = 5)

top_fits <- combined_df_clean %>%
  arrange(desc(logLik)) %>%
  head(10)

top_fits$R0 <- top_fits$beta / top_fits$gamma
top_fits$R0_lower <- top_fits$beta_lower / top_fits$gamma_lower
top_fits$R0_upper <- top_fits$beta_upper / top_fits$gamma_upper
top_fits$algorithm <- "Neldar-Mead"
top_fits$start_time <- 48
save(top_fits, file = "data/topFits_nel_s48_2021.RData")
load("data/topFits_nel_s48_2021.RData")


