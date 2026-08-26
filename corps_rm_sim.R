rm(list=ls())
library(deSolve)
library(ggplot2)
library(purrr)
library(dplyr)

mod_corps_NoRM <- function(times, outputV, parms){
  with(as.list(c(outputV, parms)), {
    arriv_rate = eval(D(expression(L / (1 + exp(-(a + k * x)))), "x"), 
                      list(L = L, a = a, k = k, x = times))
    
    dS = arriv_rate - beta1 * S * I / (S+I+R) - beta2 * S * D_inf / (S+I+R)
    
    dI = beta1 * S * I / (S+I+R) + beta2 * S * D_inf / (S+I+R) - gamma * I
    
    dR = (1-mu) * gamma * I
    
    dD_inf = mu * gamma * I - gamma_corps * D_inf
    
    #dD_deg = gamma_corps * D_inf
    
    #dN = arriv_rate - mu * gamma * I
    
    dinf_corps = beta2 * S * D_inf / (S+I+R)
    
    dinf_alive = beta1 * S * I / (S+I+R)
    
    dcum_death = mu * gamma * I
    
    #res <- c(dS, dI, dR, dD_inf, dD_deg, dN, dinf_corps, dinf_alive)
    res <- c(dS, dI, dR, dD_inf, dinf_corps, dinf_alive, dcum_death)
    
    return(list(res))
    
  })
}
mod_corps_RM <- function(times, outputV, parms){
  with(as.list(c(outputV, parms)), {
    arriv_rate = eval(D(expression(L / (1 + exp(-(a + k * x)))), "x"), 
                      list(L = L, a = a, k = k, x = times))
    
    dS = arriv_rate - beta1 * S * I / (S+I+R)
    
    dI = beta1 * S * I / (S+I+R) - gamma * I
    
    dR = (1-mu) * gamma * I
    
    dD_inf = 0
    
    #dD_deg = gamma_corps * D_inf
    
    #dN = arriv_rate - mu * gamma * I
    
    dinf_corps = 0
    
    dinf_alive = beta1 * S * I / (S+I+R)
    
    dcum_death = mu * gamma * I
    
    #res <- c(dS, dI, dR, dD_inf, dD_deg, dN, dinf_corps, dinf_alive)
    res <- c(dS, dI, dR, dD_inf, dinf_corps, dinf_alive, dcum_death)
    
    return(list(res))
    
  })
}
load("data/arrival_func_params.RData")
load("data/2022.RData")
morta22[morta22$julien_date<47, "Death_Present"] <- 0
morta22[morta22$julien_date>120, "Death_Present"] <- 0
morta22[morta22$julien_date==75, "Death_Present"] <- morta22[morta22$julien_date==75, "cumdeath"] - as.numeric(morta22[morta22$julien_date==75, "Corpse_removal"])
morta22[morta22$julien_date==76, "Death_Present"] <- morta22[morta22$julien_date==75, "Death_Present"] - as.numeric(morta22[morta22$julien_date==76, "Corpse_removal"])
morta22[morta22$julien_date==77, "Death_Present"] <- morta22[morta22$julien_date==76, "Death_Present"] - as.numeric(morta22[morta22$julien_date==77, "Corpse_removal"])
morta22[morta22$julien_date==82, "Death_Present"] <- morta22[morta22$julien_date==79, "Death_Present"] - as.numeric(morta22[morta22$julien_date==82, "Corpse_removal"])
morta22[morta22$julien_date==100, "Death_Present"] <- morta22[morta22$julien_date==95, "Death_Present"] - as.numeric(strsplit(as.character(morta22[morta22$julien_date==100, "Corpse_removal"]), split = " ", fixed = TRUE)[[1]][1])
morta22[morta22$julien_date==104, "Death_Present"] <- morta22[morta22$julien_date==100, "Death_Present"] - as.numeric(morta22[morta22$julien_date==104, "Corpse_removal"])
morta22 <- morta22 %>% filter(julien_date < 120 & julien_date > 36)

cum_death_diff <- data.frame()
nrow <- 1
# simulation with lower Ro
for(i in seq(1,7)){
  for(ratio in seq(0.2, 7, by = 0.2)){
    for(inf_period in seq(1.62, 11.62, by = 2)){
      for(inf_period_ratio in c(0.25, 0.5, 1, 1.5, 2)){
        bfitparms2 <- c(beta1 = i/inf_period, gamma = 1/inf_period, mu = 0.812, 
                        beta2 = i*ratio*inf_period_ratio/inf_period, 
                        gamma_corps = inf_period_ratio/inf_period,
                        L = logistic_func[logistic_func$year == 2022, "L"], 
                        a = logistic_func[logistic_func$year == 2022, "a"],
                        k = logistic_func[logistic_func$year == 2022, "k"])
        init_pop_size <- as.numeric(logistic_func[logistic_func$year == 2022, "L"] / (1 + exp(-(logistic_func[logistic_func$year == 2022, "a"] + logistic_func[logistic_func$year == 2022, "k"] * 46))))
        init_pop <- c(S = init_pop_size * (1 - 0.028), I = init_pop_size * 0.028, 
                      R = 0, D_inf = 0, inf_corps = 0, inf_alive = 0, cum_death = 0)
        sim_norm <- as.data.frame(rk4(init_pop, seq(46, 120), mod_corps_NoRM, bfitparms2))
        sim_rm <- as.data.frame(rk4(init_pop, seq(46, 120), mod_corps_RM, bfitparms2))
        cum_death_diff[nrow,"R_o_live"] <- i
        cum_death_diff[nrow, "R_o_corps"] <- i*ratio
        cum_death_diff[nrow, "inf_period"] <- inf_period
        cum_death_diff[nrow, "inf_period_ratio"] <- inf_period_ratio
        cum_death_diff[nrow, "death_diff"] <- tail(sim_norm,1)$cum_death - tail(sim_rm,1)$cum_death
        nrow <- nrow + 1
      }
    }
  }
}

cum_death_diff$R_o_ratio <- cum_death_diff$R_o_corps/cum_death_diff$R_o_live
cum_death_diff_plot <- cum_death_diff %>% filter(inf_period == 3.62 | inf_period == 7.62 | inf_period == 11.62)
ggplot(data = cum_death_diff %>% filter(inf_period_ratio==1), 
       aes(x = R_o_ratio, y = death_diff)) +
  facet_wrap(~inf_period, labeller = 
               labeller(inf_period = ~ paste("Infectious period: ", .x, "days"))) + 
  geom_point(aes(color = factor(R_o_live))) +
  geom_line(aes(color = factor(R_o_live))) +
  theme_bw() +
  scale_y_log10() +
  labs(y = "Difference between carcass removal vs no removal\nin total mortality number", 
       x = "Ratio of R0 (carcasses) vs R0 (live birds)",
       color = "R0 (live birds)", title = "When the infectious period of carcasses is same as that of live birds")
ggsave("../GR_DP_HPAI/figures/corps_rm_sim.png")

ggplot(data = cum_death_diff %>% filter(inf_period_ratio==0.25), 
       aes(x = R_o_ratio, y = death_diff)) +
  facet_wrap(~inf_period, labeller = 
               labeller(inf_period = ~ paste("Infectious period: ", .x, "days"))) + 
  geom_point(aes(color = factor(R_o_live))) +
  geom_line(aes(color = factor(R_o_live))) +
  theme_bw() +
  scale_y_log10() +
  labs(y = "Difference between carcass removal vs no removal\nin total mortality number", 
       x = "Ratio of R0 (carcasses) vs R0 (live birds)",
       color = "R0 (live birds)", title = "When the infectious period of carcasses is 4x longer than that of live birds")
ggsave("../GR_DP_HPAI/figures/corps_rm_sim_inf_period_4x.png")

ggplot(data = cum_death_diff_plot %>% filter(inf_period_ratio==0.25 & inf_period == 11.62), 
       aes(x = R_o_ratio, y = death_diff)) +
  #geom_point(aes(color = factor(R_o_live))) +
  geom_line(aes(color = factor(R_o_live))) +
  theme_bw() +
  scale_y_log10() +
  labs(y = "Mortality reduction by carcass removal", 
       x = "Ratio of R0 (carcasses) vs R0 (live birds)",
       color = "R0 (live birds)")
ggsave("../GR_DP_HPAI/figures/corps_rm_sim_inf_period_4x_pre.svg")

ggplot(data = cum_death_diff %>% filter(inf_period_ratio==2), 
       aes(x = R_o_ratio, y = death_diff)) +
  facet_wrap(~inf_period, labeller = 
               labeller(inf_period = ~ paste("Infectious period: ", .x, "days"))) + 
  geom_point(aes(color = factor(R_o_live))) +
  geom_line(aes(color = factor(R_o_live))) +
  theme_bw() +
  scale_y_log10() +
  labs(y = "Difference between carcass removal vs no removal\nin total mortality number", 
       x = "Ratio of R0 (carcasses) vs R0 (live birds)",
       color = "R0 (live birds)", title = "When the infectious period of carcasses is 2x shorter than that of live birds")
ggsave("../GR_DP_HPAI/figures/corps_rm_sim_inf_period_2x_shorter.png")

bfitparms2 <- c(beta1 = 2/11.62, gamma = 1/11.62, mu = 0.812, 
               beta2 = 2/11.62, gamma_corps = 1/11.62,
               L = logistic_func[logistic_func$year == 2022, "L"], 
               a = logistic_func[logistic_func$year == 2022, "a"],
               k = logistic_func[logistic_func$year == 2022, "k"])
init_pop_size <- as.numeric(logistic_func[logistic_func$year == 2022, "L"] / (1 + exp(-(logistic_func[logistic_func$year == 2022, "a"] + logistic_func[logistic_func$year == 2022, "k"] * 46))))
init_pop <- c(S = init_pop_size * (1 - 0.028), I = init_pop_size * 0.028, 
              R = 0, D_inf = 0, inf_corps = 0, inf_alive = 0, cum_death = 0)
sim_norm <- as.data.frame(rk4(init_pop, seq(46, 120), mod_corps_NoRM, bfitparms2))
sim_rm <- as.data.frame(rk4(init_pop, seq(46, 120), mod_corps_RM, bfitparms2))

# simulating realistic carcass removal when varying R_o of alive bird infection
# and keeping the carcass transmission parameters the same
START_DAY <- 75
SECOND <- 76
THIRD <- 77
FOURTH <- 79
FIFTH <- 82
SIXTH <- 100 
SEVENTH <- 104

last1 <- sim_norm %>% filter(time == START_DAY)
init_pop2 <- c(S = last1$S, I = last1$I, R = last1$R, 
                   D_inf = 0,
                   inf_corps = last1$inf_corps,
                   inf_alive = last1$inf_alive, cum_death = last1$cum_death)
sim2<-as.data.frame(rk4(init_pop2, seq(START_DAY, SECOND), mod_corps_NoRM, bfitparms2))
last2 <- tail(sim2, 1)
sim2 <- sim2 %>%
  filter(time>START_DAY)
init_pop3 <- c(S = last2$S, I = last2$I, R = last2$R, 
               # D_inf = last2$D_inf, 
               D_inf = 0,
               # D_deg = last2$D_deg, 
               N =last2$N, inf_corps = last2$inf_corps, 
               inf_alive = last2$inf_alive, cum_death = last2$cum_death)
sim3 <- as.data.frame(rk4(init_pop3, seq(SECOND, THIRD), mod_corps_NoRM, bfitparms2))
last3 <- tail(sim3, 1)
sim3 <- sim3 %>%
  filter(time>SECOND)
init_pop4 <- c(S = last3$S, I = last3$I, R = last3$R, 
               # D_inf = last3$D_inf, 
               D_inf = 0,
               # D_deg = last3$D_deg, 
               N = last3$N, inf_corps = last3$inf_corps, 
               inf_alive = last3$inf_alive,cum_death = last3$cum_death)
sim4 <- as.data.frame(rk4(init_pop4, seq(THIRD, FOURTH), mod_corps_NoRM, bfitparms2))
last4 <- tail(sim4, 1)
sim4 <- sim4 %>%
  filter(time>THIRD)
init_pop5 <- c(S = last4$S, I = last4$I, R = last4$R, 
               # D_inf = last4$D_inf, 
               D_inf = 0,
               # D_deg = last4$D_deg, 
               inf_corps = last4$inf_corps, 
               inf_alive = last4$inf_alive, 
               cum_death = last4$cum_death)
sim5 <- as.data.frame(rk4(init_pop5, seq(FOURTH, FIFTH), mod_corps_NoRM, bfitparms2))
last5 <- tail(sim5, 1)
sim5 <- sim5 %>%
  filter(time>FOURTH)

init_pop6 <- c(S = last5$S, I = last5$I, R = last5$R, 
               # D_inf = last5$D_inf, 
               D_inf = 0,
               # D_deg = last5$D_deg, 
               inf_corps = last5$inf_corps, inf_alive = last5$inf_alive,
               cum_death = last5$cum_death)
sim6 <- as.data.frame(rk4(init_pop6, seq(FIFTH, SIXTH), mod_corps_NoRM, bfitparms2))
last6 <- tail(sim6, 1)
sim6 <- sim6 %>%
  filter(time>FIFTH)

init_pop7 <- c(S = last6$S, I = last6$I, R = last6$R, 
               # D_inf = last6$D_inf,
               D_inf = 0,
               # D_deg = last6$D_deg, 
               inf_corps = last6$inf_corps, inf_alive = last6$inf_alive,
               cum_death = last6$cum_death)
sim7 <- as.data.frame(rk4(init_pop7, seq(SIXTH, SEVENTH), mod_corps_NoRM, bfitparms2))
last7 <- tail(sim7, 1)
sim7 <- sim7 %>%
  filter(time>SIXTH)

init_pop8 <- c(S = last7$S, I = last7$I, R = last7$R, 
               # D_inf = last7$D_inf, 
               D_inf = 0,
               # D_deg = last7$D_deg, 
               inf_corps = last7$inf_corps, inf_alive = last7$inf_alive,
               cum_death = last7$cum_death)
sim8 <- as.data.frame(rk4(init_pop8, seq(SEVENTH, 120), mod_corps_NoRM, bfitparms2))
sim8 <- sim8 %>%
  filter(time>SEVENTH)

sim_rm <- rbind(sim_norm %>% filter(time < START_DAY), sim2, sim3, sim4, sim5, sim6, sim7, sim8)

# assuming no corpse transmission occuring after removing carcasses
sim2_rm <-as.data.frame(rk4(init_pop2, seq(START_DAY, FIFTH), mod_corps_RM, bfitparms2))
last2_rm <- tail(sim2_rm, 1)

init_pop3_rm <- c(S = last2_rm$S, I = last2_rm$I, R = last2_rm$R, 
               D_inf = last2_rm$D_inf, inf_corps = last2_rm$inf_corps, 
               inf_alive = last2_rm$inf_alive, cum_death = last2_rm$cum_death)
sim3_rm <- as.data.frame(rk4(init_pop3_rm, seq(FIFTH, SIXTH), mod_corps_NoRM, bfitparms2))
last3_rm <- tail(sim3_rm, 1)

init_pop4_rm <- c(S = last3_rm$S, I = last3_rm$I, R = last3_rm$R, 
                  D_inf = 0, inf_corps = last3_rm$inf_corps, 
                  inf_alive = last3_rm$inf_alive, cum_death = last3_rm$cum_death)
sim4_rm <- as.data.frame(rk4(init_pop4_rm, seq(SIXTH, SEVENTH), mod_corps_RM, bfitparms2))
last4_rm <- tail(sim4_rm, 1)

init_pop8_rm <- c(S = last4_rm$S, I = last4_rm$I, R = last4_rm$R, 
               # D_inf = last7$D_inf, 
               D_inf = 0,
               # D_deg = last7$D_deg, 
               inf_corps = last4_rm$inf_corps, inf_alive = last4_rm$inf_alive,
               cum_death = last4_rm$cum_death)
sim8_rm <- as.data.frame(rk4(init_pop8_rm, seq(SEVENTH, 120), mod_corps_NoRM, bfitparms2))
sim8_rm <- sim8_rm %>%
  filter(time>SEVENTH)

sim_rm <- rbind(sim2_rm, sim3_rm, sim4_rm, sim8_rm)
