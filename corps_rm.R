rm(list=ls())
library(deSolve)
library(ggplot2)
library(purrr)
library(dplyr)

mod_corps_NoRM <- function(times, outputV, parms){
  with(as.list(c(outputV, parms)), {
    arriv_rate = eval(D(expression(L / (1 + exp(-(a + k * x)))), "x"), 
                      list(L = L, a = a, k = k, x = times))
    
    dS = arriv_rate - beta1 * S * I / N - beta2 * S * D_inf / N
    
    dI = beta1 * S * I / N + beta2 * S * D_inf / N - gamma * I
    
    dR = (1-mu) * gamma * I
    
    dD_inf = mu * gamma * I - gamma_corps * D_inf
    
    dD_deg = gamma_corps * D_inf
    
    dN = arriv_rate - mu * gamma * I
    
    dinf_corps = beta2 * S * D_inf / N
    
    dinf_alive = beta1 * S * I / N
    
    res <- c(dS, dI, dR, dD_inf, dD_deg, dN, dinf_corps, dinf_alive)
    
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

param.grid <- expand.grid(
  beta1_2_ratio = seq(0, 1, by = 0.01),
  gamma_corps_inv = seq(1, 40, by = 1)
)



START_DAY <- 75
SECOND <- 76
THIRD <- 77
FOURTH <- 79
FIFTH <- 82
SIXTH <- 100
SEVENTH <- 104

# if moving 25 days earlier
START_DAY <- START_DAY - 25
SECOND <- SECOND - 25
THIRD <- THIRD - 25
FOURTH <- FOURTH - 25
FIFTH <- FIFTH - 25
SIXTH <- SIXTH - 25
SEVENTH <- SEVENTH - 25

# fit the model to the first time period without corpse removal
simdf_ts <- map(1:nrow(param.grid), function(idx) {
  parms22 <- c(beta1 = 7.07/11.62, gamma = 1 / 11.62, mu = 0.812, 
               beta2 = 7.07/11.62 * param.grid$beta1_2_ratio[idx],
               gamma_corps = 1/ param.grid$gamma_corps_inv[idx],
               L = logistic_func[logistic_func$year == 2022, "L"], 
               a = logistic_func[logistic_func$year == 2022, "a"],
               k = logistic_func[logistic_func$year == 2022, "k"])
  init_pop_size <- as.numeric(parms22["L"] / (1 + exp(-(parms22["a"] + parms22["k"] * 46))))
  init_pop <- c(S = init_pop_size * (1 - 0.028),
                I = init_pop_size * 0.028, R = 0, D_inf = 0, D_deg = 0, N = init_pop_size, inf_corps = 0, inf_alive = 0)
  sim <- as.data.frame(rk4(init_pop, seq(46, 75), mod_corps_NoRM, parms22))
  sim_sub <- sim %>% filter(time %in% morta22$julien_date)
  morta22_sub <- morta22 %>% filter(julien_date %in% sim_sub$time)
  mse_cum <- sum((sim_sub$D_inf + sim_sub$D_deg - morta22_sub$cumdeath)^2) / sum((morta22_sub$cumdeath)^2)
  sim$mse_cum <- mse_cum
  print(paste0(mse_cum))
  return(sim)
})

simdf1 <- bind_rows(simdf_ts, .id = "id")
idx_min <- simdf1$id[which.min(simdf1$mse_cum)]

bfitparms <- c(beta1 = 7.07/11.62, gamma = 1 / 11.62, mu = 0.812, beta2 = 7.07/11.62 * param.grid$beta1_2_ratio[as.numeric(idx_min)],
               gamma_corps = 1/ param.grid$gamma_corps_inv[as.numeric(idx_min)],
               L = logistic_func[logistic_func$year == 2022, "L"], a = logistic_func[logistic_func$year == 2022, "a"],
               k = logistic_func[logistic_func$year == 2022, "k"])
simdf1_bfit <- simdf1 %>% filter(id == idx_min)

last1 <- simdf1_bfit %>% filter(time == START_DAY)
# simulate later stage of corpse removal using the best fitted parameters of corpse transmission
last1$D <- last1$D_inf + last1$D_deg
# last1$N <- last1$N - 202 * last1$D_deg/last1$D
last1$D_inf <- last1$D_inf - 202 * last1$D_inf/last1$D
last1$D_deg <- last1$D_deg - 202 * last1$D_deg/last1$D
sim1 <- simdf1_bfit %>% select(-id, -mse_cum)
sim1$D <- sim1$D_deg + sim1$D_inf

init_pop2 <- c(S = last1$S, I = last1$I, R = last1$R, D_inf = last1$D_inf, 
               D_deg = last1$D_deg, N = last1$N, inf_corps = last1$inf_corps, inf_alive = last1$inf_alive)
sim2<-as.data.frame(rk4(init_pop2, seq(START_DAY, SECOND), mod_corps_NoRM, bfitparms))
last2 <- tail(sim2, 1)
# sim2$D_inf <- sim2$D_inf + 202 * last1$D_inf/last1$D
# sim2$D_deg <- sim2$D_deg + 202 * last1$D_deg/last1$D
last2$D <- last2$D_inf + last2$D_deg
# last2$N <- last2$N - 243 * last2$D_deg/last2$D
if((last2$D_inf - 243 * last2$D_inf/last2$D) > 0){
  last2$D_inf <- last2$D_inf - 243 * last2$D_inf/last2$D  
} else {
  last2$D_inf <- 0
}
if((last2$D_deg - 243 * last2$D_deg/last2$D) > 0){
  last2$D_deg <- last2$D_deg - 243 * last2$D_deg/last2$D
} else {
  last2$D_deg <- 0
}
sim2 <- sim2 %>%
  filter(time>START_DAY)
sim2$D <- sim2$D_deg + sim2$D_inf + 202

init_pop3 <- c(S = last2$S, I = last2$I, R = last2$R, D_inf = last2$D_inf, 
               D_deg = last2$D_deg, N =last2$N, inf_corps = last2$inf_corps, inf_alive = last2$inf_alive)
sim3 <- as.data.frame(rk4(init_pop3, seq(SECOND, THIRD), mod_corps_NoRM, bfitparms))
last3 <- tail(sim3, 1)
# sim3$D_inf <- sim3$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D
# sim3$D_deg <- sim3$D_deg + 202 * last1$D_deg/last1$D + 243 * last2$D_deg/last2$D
last3$D <- last3$D_inf + last3$D_deg
# last3$N <- last3$N - 300 * last3$D_deg/last3$D
if((last3$D_inf - 300 * last3$D_inf/last3$D) > 0){
  last3$D_inf <- last3$D_inf - 300 * last3$D_inf/last3$D 
} else {
  last3$D_inf <- 0
}
if((last3$D_deg - 300 * last3$D_deg/last3$D) > 0){
  last3$D_deg <- last3$D_deg - 300 * last3$D_deg/last3$D  
} else {
  last3$D_deg <- 0
}
sim3 <- sim3 %>%
  filter(time>SECOND)
sim3$D <- sim3$D_deg + sim3$D_inf + 202 + 243

init_pop4 <- c(S = last3$S, I = last3$I, R = last3$R, D_inf = last3$D_inf, 
               D_deg = last3$D_deg, N = last3$N, inf_corps = last3$inf_corps, inf_alive = last3$inf_alive)
sim4 <- as.data.frame(rk4(init_pop4, seq(THIRD, FOURTH), mod_corps_NoRM, bfitparms))
last4 <- tail(sim4, 1)
# sim4$D_inf <- sim4$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D
# sim4$D_deg <- sim4$D_deg + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_deg/last3$D
last4$D <- last4$D_inf + last4$D_deg
# last4$N <- last4$N - 226 * last4$D_deg/last4$D
if((last4$D_inf - 226 * last4$D_inf/last4$D) > 0){
  last4$D_inf <- last4$D_inf - 226 * last4$D_inf/last4$D  
} else {
  last4$D_inf <- 0
}
if((last4$D_deg - 226 * last4$D_deg/last4$D) > 0) {
  last4$D_deg <- last4$D_deg - 226 * last4$D_deg/last4$D  
} else {
  last4$D_deg <- 0
}
sim4 <- sim4 %>%
  filter(time>THIRD)
sim4$D <- sim4$D_deg + sim4$D_inf + 202 + 243 + 300

init_pop5 <- c(S = last4$S, I = last4$I, R = last4$R, D_inf = last4$D_inf, 
               D_deg = last4$D_deg, N = last4$N, inf_corps = last4$inf_corps, inf_alive = last4$inf_alive)
sim5 <- as.data.frame(rk4(init_pop5, seq(FOURTH, FIFTH), mod_corps_NoRM, bfitparms))
last5 <- tail(sim5, 1)
# sim5$D_inf <- sim5$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D
# sim5$D_deg <- sim5$D_deg + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D
last5$D <- last5$D_inf + last5$D_deg
# last5$N <- last5$N - 172 * last5$D_deg/last5$D
if((last5$D_inf - 172 * last5$D_inf/last5$D)>0){
  last5$D_inf <- last5$D_inf - 172 * last5$D_inf/last5$D  
}else{
  last5$D_inf <- 0
}
if((last5$D_deg - 172 * last5$D_deg/last5$D) > 0){
  last5$D_deg <- last5$D_deg - 172 * last5$D_deg/last5$D  
} else {
  last5$D_deg <- 0
}
sim5 <- sim5 %>%
  filter(time>FOURTH)
sim5$D <- sim5$D_deg + sim5$D_inf + 202 + 243 + 300 + 226

init_pop6 <- c(S = last5$S, I = last5$I, R = last5$R, D_inf = last5$D_inf, D_deg = last5$D_deg, N = last5$N,
               inf_corps = last5$inf_corps, inf_alive = last5$inf_alive)
sim6 <- as.data.frame(rk4(init_pop6, seq(FIFTH, SIXTH), mod_corps_NoRM, bfitparms))
last6 <- tail(sim6, 1)
# sim6$D_inf <- sim6$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D
# sim6$D_deg <- sim6$D_deg + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D
last6$D <- last6$D_inf + last6$D_deg
# last6$N <- last6$N - 110 * last6$D_deg/last6$D
if((last6$D_inf - 110 * last6$D_inf/last6$D) > 0){
  last6$D_inf <- last6$D_inf - 110 * last6$D_inf/last6$D
} else {
  last6$D_inf <- 0
}
if((last6$D_deg - 110 * last6$D_deg/last6$D) >0){
  last6$D_deg <- last6$D_deg - 110 * last6$D_deg/last6$D
} else {
  last6$D_deg <- 0
}
sim6 <- sim6 %>%
  filter(time>FIFTH)
sim6$D <- sim6$D_deg + sim6$D_inf + 202 + 243 + 300 + 226 + 172

init_pop7 <- c(S = last6$S, I = last6$I, R = last6$R, D_inf = last6$D_inf, D_deg = last6$D_deg, N = last6$N,
               inf_corps = last6$inf_corps, inf_alive = last6$inf_alive)
sim7 <- as.data.frame(rk4(init_pop7, seq(SIXTH, SEVENTH), mod_corps_NoRM, bfitparms))
last7 <- tail(sim7, 1)
# sim7$D_inf <- sim7$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D + 110 * last6$D_inf/last6$D
# sim7$D_deg <- sim7$D_deg + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D + 110 * last6$D_inf/last6$D
last7$D <- last7$D_inf + last7$D_deg
# last7$N <- last7$N - 167 * last7$D_deg/last7$D
if((last7$D_inf - 167 * last7$D_inf/last7$D) > 0){
  last7$D_inf <- last7$D_inf - 167 * last7$D_inf/last7$D
} else {
  last7$D_inf <- 0
}
if((last7$D_deg - 167 * last7$D_deg/last7$D) >0){
  last7$D_deg <- last7$D_deg - 167 * last7$D_deg/last7$D
} else {
  last7$D_deg <- 0
}
sim7 <- sim7 %>%
  filter(time>SIXTH)
sim7$D <- sim7$D_deg + sim7$D_inf + 202 + 243 + 300 + 226 + 172 + 110

init_pop8 <- c(S = last7$S, I = last7$I, R = last7$R, D_inf = last7$D_inf, 
               D_deg = last7$D_deg, N = last7$N, inf_corps = last7$inf_corps, 
               inf_alive = last7$inf_alive)
sim8 <- as.data.frame(rk4(init_pop8, seq(SEVENTH, 120), mod_corps_NoRM, bfitparms))
# sim8$D_inf <- sim8$D_inf + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D + 110 * last6$D_inf/last6$D + 167 * last7$D_inf/last7$D
# sim8$D_deg <- sim8$D_deg + 202 * last1$D_inf/last1$D + 243 * last2$D_inf/last2$D + 300 * last3$D_inf/last3$D + 226 * last4$D_inf/last4$D + 172 * last5$D_inf/last5$D + 110 * last6$D_inf/last6$D + 167 * last7$D_inf/last7$D
sim8 <- sim8 %>%
  filter(time>SEVENTH)
sim8$D <- sim8$D_deg + sim8$D_inf + 202 + 243 + 300 + 226 + 172 + 110 + 167

sim_rm <- rbind(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8)
init_pop_size <- as.numeric(logistic_func[logistic_func$year == 2022, "L"] / (1 + exp(-(logistic_func[logistic_func$year == 2022, "a"] + logistic_func[logistic_func$year == 2022, "k"] * 46))))
init_pop <- c(S = init_pop_size * (1 - 0.028),
              I = init_pop_size * 0.028, R = 0, D_inf = 0, D_deg = 0, N = init_pop_size, inf_corps = 0, inf_alive = 0)
# init_pop_size <- as.numeric(logistic_func[logistic_func$year == 2021, "L"] / (1 + exp(-(logistic_func[logistic_func$year == 2021, "a"] + logistic_func[logistic_func$year == 2021, "k"] * 48))))
# init_pop <- c(S = init_pop_size * (1 - 0.0003),
#               I = init_pop_size * 0.0003, R = 0, D_inf = 0, D_deg = 0, N = init_pop_size, inf_corps = 0, inf_alive = 0)
sim_norm <- as.data.frame(rk4(init_pop, seq(46, 120), mod_corps_NoRM, bfitparms))

ggplot(data = sim_rm, aes(x = time)) +
  geom_line(aes(y = S, color = "Susceptible")) +
  geom_line(aes(y = inf_corps, color = "Cumulative infections from dead birds", linetype = "with corpse removal")) +
  geom_line(aes(y = inf_alive, color = "Cumulative infections from alive birds", linetype = "with corpse removal")) +
  # geom_line(aes(y = I, color = "Infected")) +
  geom_line(aes(y = D, color = "Cumulative dead birds", linetype = "with corpse removal")) +
  labs(x = "Time (julian date)", y = "Population size", color = "Population", fill = "95% confidence interval") +
  theme_minimal() +
  # geom_line(data = morta22, aes(x = julien_date, y = cumdeath, linetype = "Observed death")) +
  geom_line(data = sim_norm, aes(x = time, y = D_inf + D_deg, color = "Cumulative dead birds", linetype = "without corpse removal")) +
  geom_line(data = sim_norm, aes(x = time, y = inf_alive, color = "Cumulative infections from alive birds", linetype = "without corpse removal")) +
  geom_line(data = sim_norm, aes(x = time, y = inf_corps, color = "Cumulative infections from dead birds", linetype = "without corpse removal")) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_color_manual(values = c("black", "red", "darkred", "darkgrey")) +
  geom_vline(aes(xintercept = START_DAY), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = SECOND), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = THIRD), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = FOURTH), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = FIFTH), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = SIXTH), linetype = 3, color = "grey") +
  geom_vline(aes(xintercept = SEVENTH), linetype = 3, color = "grey") +
  theme(legend.title = element_blank(), legend.position = "inside", legend.position.inside = c(0.75, 0.6))
ggsave("figures/corpse_removal_impacts.png", width = 5.5, height = 5)

# Print summary of impact
cat("\n=== IMPACT ASSESSMENT ===\n")
final_deaths_with <- tail(sim_rm$D, 1)
final_deaths_without <- tail(sim_norm$D_inf + sim_norm$D_deg, 1)
deaths_prevented <- final_deaths_without - final_deaths_with

cat("Final deaths WITH corpse removal:", round(final_deaths_with, 0), "\n")
cat("Final deaths WITHOUT corpse removal:", round(final_deaths_without, 0), "\n")
cat("Deaths prevented by corpse removal:", round(deaths_prevented, 0), "\n")
cat("Percentage reduction:", round(100 * deaths_prevented / final_deaths_without, 1), "%\n\n")

# Impact on infection sources
final_inf_corps_with <- tail(sim_rm$inf_corps, 1)
final_inf_corps_without <- tail(sim_norm$inf_corps, 1)
final_inf_alive_with <- tail(sim_rm$inf_alive, 1)
final_inf_alive_without <- tail(sim_norm$inf_alive, 1)

cat("Infections from corpses WITH removal:", round(final_inf_corps_with, 0), "\n")
cat("Infections from corpses WITHOUT removal:", round(final_inf_corps_without, 0), "\n")
cat("Corpse infections prevented:", round(final_inf_corps_without - final_inf_corps_with, 0), "\n\n")

cat("Infections from live birds WITH removal:", round(final_inf_alive_with, 0), "\n")
cat("Infections from live birds WITHOUT removal:", round(final_inf_alive_without, 0), "\n")
cat("Change in live bird infections:", round(final_inf_alive_with - final_inf_alive_without, 0), "\n")

# Total corpses removed
total_removed <- sum(c(202, 243, 300, 226, 172, 110, 167))
cat("\nTotal corpses removed:", total_removed, "\n")
cat("Deaths prevented per corpse removed:", round(deaths_prevented / total_removed, 2), "\n")
