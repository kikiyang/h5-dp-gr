rm(list=ls())
library(dplyr)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(cowplot)

# simple transmission model
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

load("data/top_fits_nel_s46_2022.RData")
bestfit22 <- head(top_fits_nel_s46, 1)
load("data/topFits_nel_s48_2021.RData")
bestfit_21 <- head(top_fits, 1)
load("data/arrival_func_params.RData")
load("data/2021.RData")
params21 <- c(beta = bestfit_21$beta, gamma = bestfit_21$gamma, 
              mu = bestfit_21$mu, 
              L = logistic_func[logistic_func$year == 2021, "L"], 
              a = logistic_func[logistic_func$year == 2021, "a"], 
              k = logistic_func[logistic_func$year == 2021, "k"])
time_vec21 <- seq(48, 150)
init_pop_size21 <- as.numeric(params21["L"] / (1 + exp(-(params21["a"] + params21["k"] * min(time_vec21)))))
init_pop21 <- c(S = init_pop_size21 * (1 - bestfit_21$I0),
                I = init_pop_size21 * bestfit_21$I0, R = 0, D = 0)
sim21 <- as.data.frame(rk4(init_pop21, time_vec21, mod4, params21))

params21_lower <- c(beta = bestfit_21$beta_lower, gamma = bestfit_21$gamma_lower,
                    mu = bestfit_21$mu_lower, 
                    L = logistic_func[logistic_func$year == 2021, "L"], 
                    a = logistic_func[logistic_func$year == 2021, "a"], 
                    k = logistic_func[logistic_func$year == 2021, "k"])
init_pop21_lower <- c(S = init_pop_size21 * (1 - bestfit_21$I0_lower), 
                      I = init_pop_size21 * bestfit_21$I0_lower, R = 0, D = 0)
sim21_lower <- as.data.frame(rk4(init_pop21_lower, time_vec21, mod4, params21_lower))

params21_upper <- c(beta = bestfit_21$beta_upper, gamma = bestfit_21$gamma_upper,
                    mu = bestfit_21$mu_upper, 
                    L = logistic_func[logistic_func$year == 2021, "L"], 
                    a = logistic_func[logistic_func$year == 2021, "a"], 
                    k = logistic_func[logistic_func$year == 2021, "k"])
init_pop21_upper <- c(S = init_pop_size21 * (1 - bestfit_21$I0_upper),
                      I = init_pop_size21 * bestfit_21$I0_upper, R = 0, D = 0)
sim21_upper <- as.data.frame(rk4(init_pop21_upper, time_vec21, mod4, params21_upper))

p21 <- ggplot(data = sim21, aes(x = time)) +
  geom_ribbon(aes(ymin = sim21_lower$S, ymax = sim21_upper$S, fill = "Susceptible"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim21_lower$I, ymax = sim21_upper$I, fill = "Infected"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim21_lower$D, ymax = sim21_upper$D, fill = "Dead"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim21_lower$pop_size, ymax = sim21_upper$pop_size, fill = "Alive"), alpha = 0.2) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = D, color = "Dead")) +
  geom_line(aes(y = pop_size, color = "Alive")) +
  labs(x = "Julian date", y = "Population size", color = "Population", fill = "95% confidence interval") +
  theme_minimal() +
  geom_line(data = morta21, aes(x = julien_date, y = cum_death, color = "Observed death"), size = 1) +
  geom_line(data = nest21, aes(x = julien_date, y = est_alive, color = "Observed alive"), size = 1) +
  scale_color_manual(values = c("black", "darkred", "steelblue", "darkgrey", "pink", "darkgreen")) +
  scale_fill_manual(values = c("black", "darkred", "steelblue", "darkgreen")) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 2400))

time_vec <- seq(46, 120)
params <- c(beta = bestfit22$beta, gamma = bestfit22$gamma, 
            mu = bestfit22$mu, 
            L = logistic_func[logistic_func$year == 2022, "L"], 
            a = logistic_func[logistic_func$year == 2022, "a"], 
            k = logistic_func[logistic_func$year == 2022, "k"])
init_pop_size <- as.numeric(params["L"] / (1 + exp(-(params["a"] + params["k"] * min(time_vec)))))
init_pop <- c(S = init_pop_size * (1 - bestfit22$I0),
              I = init_pop_size * bestfit22$I0, R = 0, D = 0)
sim <- as.data.frame(rk4(init_pop, time_vec, mod4, params))

params_lower <- c(beta = bestfit22$beta_lower, gamma = bestfit22$gamma_lower,
                  mu = bestfit22$mu_lower, 
                  L = logistic_func[logistic_func$year == 2022, "L"], 
                  a = logistic_func[logistic_func$year == 2022, "a"], 
                  k = logistic_func[logistic_func$year == 2022, "k"])
init_pop_lower <- c(S = init_pop_size * (1 - bestfit22$I0_lower),
                    I = init_pop_size * bestfit22$I0_lower, R = 0, D = 0)
sim_lower <- as.data.frame(rk4(init_pop_lower, time_vec, mod4, params_lower))

params_upper <- c(beta = bestfit22$beta_upper, gamma = bestfit22$gamma_upper,
                  mu = bestfit22$mu_upper, 
                  L = logistic_func[logistic_func$year == 2022, "L"], 
                  a = logistic_func[logistic_func$year == 2022, "a"], 
                  k = logistic_func[logistic_func$year == 2022, "k"])
init_pop_upper <- c(S = init_pop_size * (1 - bestfit22$I0_upper),
                    I = init_pop_size * bestfit22$I0_upper, R = 0, D = 0)
sim_upper <- as.data.frame(rk4(init_pop_upper, time_vec, mod4, params_upper))

load("data/2022.RData")
morta22[morta22$julien_date<47, "Death_Present"] <- 0
morta22[morta22$julien_date>120, "Death_Present"] <- 0
morta22[morta22$julien_date==75, "Death_Present"] <- morta22[morta22$julien_date==75, "cumdeath"] - as.numeric(morta22[morta22$julien_date==75, "Corpse_removal"])
morta22[morta22$julien_date==76, "Death_Present"] <- morta22[morta22$julien_date==75, "Death_Present"] - as.numeric(morta22[morta22$julien_date==76, "Corpse_removal"])
morta22[morta22$julien_date==77, "Death_Present"] <- morta22[morta22$julien_date==76, "Death_Present"] - as.numeric(morta22[morta22$julien_date==77, "Corpse_removal"])
morta22[morta22$julien_date==82, "Death_Present"] <- morta22[morta22$julien_date==79, "Death_Present"] - as.numeric(morta22[morta22$julien_date==82, "Corpse_removal"])
morta22[morta22$julien_date==100, "Death_Present"] <- morta22[morta22$julien_date==95, "Death_Present"] - as.numeric(strsplit(as.character(morta22[morta22$julien_date==100, "Corpse_removal"]), split = " ", fixed = TRUE)[[1]][1])
morta22[morta22$julien_date==104, "Death_Present"] <- morta22[morta22$julien_date==100, "Death_Present"] - as.numeric(morta22[morta22$julien_date==104, "Corpse_removal"])
morta22 <- morta22 %>% filter(julien_date < 120)

p22 <- ggplot(data = sim, aes(x = time)) +
  geom_ribbon(aes(ymin = sim_lower$S, ymax = sim_upper$S, fill = "Susceptible"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim_lower$I, ymax = sim_upper$I, fill = "Infected"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim_lower$D, ymax = sim_upper$D, fill = "Dead"), alpha = 0.2) +
  geom_ribbon(aes(ymin = sim_lower$pop_size, ymax = sim_upper$pop_size, fill = "Alive"), alpha = 0.2) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = D, color = "Dead")) +
  geom_line(aes(y = pop_size, color = "Alive")) +
  labs(x = "Julian date", y = "Population size", color = "Population", fill = "95% confidence interval") +
  theme_minimal() +
  geom_line(data = morta22, aes(x = julien_date, y = cumdeath, color = "Observed death"), size = 1) +
  geom_line(data = morta22, aes(x = julien_date, y = TOTAL_DP_ALIVE, color = "Observed alive"), size = 1) +
  scale_color_manual(values = c("black", "darkred", "steelblue", "darkgrey", "pink", "darkgreen")) +
  scale_fill_manual(values = c("black", "darkred", "steelblue", "darkgreen")) +
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 2400))

combined_plot <- p21+p22 + plot_annotation(tag_levels = 'A')
ggsave(plot = combined_plot, "figures/outbreak_mle_fit_v.png", width = 10, height = 4)


