rm(list=ls())
library(deSolve)
library(dplyr)
library(ggplot2)
library(purrr)

# transmission model
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

# constant arrival rate as comparison
mod5 <- function(times, outputV, parms){
  with(as.list(c(outputV, parms)),{
    
    N = S + I + R
    
    dS = arriv_rate - beta * S * I / N
    
    dI = beta * S * I / N - gamma * I
    
    dR = (1-mu) * gamma * I
    
    dD = mu * gamma * I
    
    res <- c(dS, dI, dR, dD)
    
    return(list(res, pop_size = N))
  })
}


# arrival function
load("data/arrival_func_params.RData")
load("data/arrival_combined.RData")

# simulate the relationship between starting time and death (at end point t) with different R0 for 2021 arrivals
param.grid <- expand.grid(
  I0 = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
  R0 = seq(1, 10, by = 0.5),
  start_time = seq(0, 120, by = 5)
)

# control group: constant arrival rate
simdf_const22 <- map(1:nrow(param.grid), function(idx){
  params <- c(beta = param.grid$R0[idx] * 1/15, gamma = 1/15, mu = 0.9, arriv_rate = 1937/120)
  time_vec <- seq(param.grid$start_time[idx], 120)
  init_pop <- c(S = as.numeric(params["arriv_rate"])*param.grid$start_time[idx]*(1-param.grid$I0[idx]),
                I = as.numeric(params["arriv_rate"])*param.grid$start_time[idx]*param.grid$I0[idx],
                R = 0, D = 0)
  sim <- as.data.frame(rk4(init_pop, time_vec, mod5, params))
  return(sim)
})

sims_const22 <- bind_rows(simdf_const22, .id = "id")
sims_const22$I0 <- param.grid$I0[as.numeric(sims_const22$id)]
sims_const22$R0 <- param.grid$R0[as.numeric(sims_const22$id)]
sims_const22$start_time <- param.grid$start_time[as.numeric(sims_const22$id)]

sims_const22_episize <- sims_const22 %>% filter(time == 120)
ggplot(data = sims_const22_episize, aes(x = start_time, y = I0, z = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_contour_filled(breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 1850, 1900, 1925, 1937, 1943, 1950, 1953, 1956, 1959, 1962)) +
  labs(x = "Outbreak start day (julian date)", y = "Initial ratio of infectious individuals", 
       fill = "Final epidemic size\n(cumulative deaths)\nat day 120")

ggplot(data = sims_const22_episize %>% filter(R0<5.5), aes(x = start_time, y = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_line(aes(group = I0, color = I0)) +
  scale_color_viridis_c(name = "Initial ratio of\ninfectious\nindividuals", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/const_arriv_rate_2022.png")
ggplot(data = sims_const22_episize %>% filter(R0<5.5 & I0 %in% c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)), aes(x = start_time, y = D)) +
  facet_wrap(~I0, labeller = labeller(I0 = label_both)) +
  geom_line(aes(group = R0, color = R0)) +
  scale_color_viridis_c(name = "R0", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")

# control group 2021: constant arrival rate
simdf_const21 <- map(1:nrow(param.grid), function(idx){
  params <- c(beta = param.grid$R0[idx] * 1/15, gamma = 1/15, mu = 0.9, arriv_rate = 2237/120)
  time_vec <- seq(param.grid$start_time[idx], 120)
  init_pop <- c(S = as.numeric(params["arriv_rate"])*param.grid$start_time[idx]*(1-param.grid$I0[idx]),
                I = as.numeric(params["arriv_rate"])*param.grid$start_time[idx]*param.grid$I0[idx],
                R = 0, D = 0)
  sim <- as.data.frame(rk4(init_pop, time_vec, mod5, params))
  return(sim)
})

sims_const21 <- bind_rows(simdf_const21, .id = "id")
sims_const21$I0 <- param.grid$I0[as.numeric(sims_const21$id)]
sims_const21$R0 <- param.grid$R0[as.numeric(sims_const21$id)]
sims_const21$start_time <- param.grid$start_time[as.numeric(sims_const21$id)]

sims_const21_episize <- sims_const21 %>% filter(time == 120)
ggplot(data = sims_const21_episize, aes(x = start_time, y = I0, z = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_contour_filled(breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1450, 1500, 1550, 1600, 1625, 1650, 1675, 1700, 1703, 1706, 1709, 1712)) +
  labs(x = "Outbreak start day (julian date)", y = "Initial ratio of infectious individuals", 
       fill = "Final epidemic size\n(cumulative deaths)\nat day 120")

ggplot(data = sims_const21_episize %>% filter(R0<5.5), aes(x = start_time, y = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_line(aes(group = I0, color = I0)) +
  scale_color_viridis_c(name = "Initial ratio of\ninfectious\nindividuals", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/constant_arriv_rate_2021.png")
ggplot(data = sims_const21_episize %>% filter(I0 %in% c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)), aes(x = start_time, y = D)) +
  facet_wrap(~I0, labeller = labeller(I0 = label_both)) +
  geom_line(aes(group = R0, color = R0)) +
  scale_color_viridis_c(name = "R0", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")

const21_epi_day50 <- sims_const21_episize %>%
  filter(start_time==50) %>%
  group_by(R0, I0)

const21_dif <- sims_const21_episize %>%
  left_join(const21_epi_day50, by = c("R0", "I0")) %>%
  mutate(D_dif = D.y - D.x) %>%
  select(I0, R0, D_dif, D.x, D.y, start_time.x)

ggplot(const21_dif, aes(x = start_time.x, y = I0, z = D_dif)) +
  facet_wrap(~R0) +
  geom_contour_filled()

## arrivals as logistic growth
## 2021
simdf_mod_ts21 <- map(1:nrow(param.grid), function(idx) {
  params <- c(beta = param.grid$R0[idx] * 1/15, gamma = 1/15, mu = 0.9,
              L = logistic_func[logistic_func$year == 2021, "L"], 
              a = logistic_func[logistic_func$year == 2021, "a"], 
              k = logistic_func[logistic_func$year == 2021, "k"])
  time_vec <- seq(param.grid$start_time[idx], 120)
  # based on the start time, the susceptible number is different
  init_pop_size <- as.numeric(params["L"] / (1 + exp(-(params["a"] + params["k"] * param.grid$start_time[idx]))))
  init_pop <- c(S = init_pop_size * (1 - param.grid$I0[idx]), 
                I = init_pop_size * param.grid$I0[idx], R = 0, D = 0)
  sim <- as.data.frame(rk4(init_pop, time_vec, mod4, params))
  return(sim)
})

simdfs21 <- bind_rows(simdf_mod_ts21, .id = "id")
simdfs21$I0 <- param.grid$I0[as.numeric(simdfs21$id)]
simdfs21$R0 <- param.grid$R0[as.numeric(simdfs21$id)]
simdfs21$start_time <- param.grid$start_time[as.numeric(simdfs21$id)]

simdfs_epi_size21 <- simdfs21 %>% filter(time == 120)
ggplot(data = simdfs_epi_size21 %>% filter(R0<5.5), aes(x = start_time, y = I0, z = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_contour_filled(breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 1850, 1900, 1925, 1937, 1943, 1950, 1953, 1956, 1959, 1962)) +
  labs(x = "Outbreak start day (julian date)", y = "Initial ratio of infectious individuals", 
       fill = "Final epidemic size\n(cumulative deaths)\nat day 120")
ggsave("figures/2021arriv_start_I0.png")

simdfs_epi_size_I1_21 <- simdfs_epi_size21 %>% filter(time == 120)
ggplot(data = simdfs_epi_size_I1_21 %>% filter(R0<5.5), aes(x = start_time, y = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_line(aes(group = I0, color = I0)) +
  scale_color_viridis_c(name = "Initial ratio of\ninfectious\nindividuals", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/2021_start_D_R0_facets.png", width = 7, height = 5)

ggplot(data = simdfs_epi_size_I1_21 %>% filter(R0<5.5 & I0 %in% c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)), aes(x = start_time, y = D)) +
  facet_wrap(~I0, labeller = labeller(I0 = label_both)) +
  geom_line(aes(group = R0, color = R0)) +
  scale_color_viridis_c(name = "R0", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/2021_start_D_I0_facets.png", width = 7, height = 5)

# Find the start_time that gives maximum D for each R0 and I0 combination
max_death_summary21 <- simdfs_epi_size_I1_21 %>%
  group_by(R0, I0) %>%
  summarise(
    max_death = max(D),
    optimal_start_time = start_time[which.max(D)],
    .groups = 'drop'
  ) %>%
  arrange(R0, I0)

# View the summary
print(max_death_summary21)

# Create a heatmap showing optimal start times
ggplot(max_death_summary21, aes(x = R0, y = I0, fill = optimal_start_time)) +
  geom_tile() +
    scale_fill_gradient2(name = "Optimal\nStart Time", 
                       low = "blue", mid = "white", high = "red", 
                       midpoint = 50) +
  labs(x = "R0", y = "Initial Infectious Proportion (I0)", 
       title = "Optimal Start Time for Maximum Deaths") +
  theme_minimal()

# Create a heatmap showing maximum deaths achieved
ggplot(max_death_summary21, aes(x = R0, y = I0, fill = max_death)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Maximum\nDeaths") +
  labs(x = "R0", y = "Initial Infectious Proportion (I0)", 
       title = "Maximum Deaths Achieved") +
  theme_minimal()

## arrivals as logistic growth
## 2022
simdf_mod_ts22 <- map(1:nrow(param.grid), function(idx) {
  params <- c(beta = param.grid$R0[idx] * 1/15, gamma = 1/15, mu = 0.9,
              L = logistic_func[logistic_func$year == 2022, "L"], 
              a = logistic_func[logistic_func$year == 2022, "a"], 
              k = logistic_func[logistic_func$year == 2022, "k"])
  time_vec <- seq(param.grid$start_time[idx], 120)
  init_pop_size <- as.numeric(params["L"] / (1 + exp(-(params["a"] + params["k"] * param.grid$start_time[idx]))))
  init_pop <- c(S = init_pop_size * (1 - param.grid$I0[idx]), 
                I = init_pop_size * param.grid$I0[idx], R = 0, D = 0)
  sim <- as.data.frame(rk4(init_pop, time_vec, mod4, params))
  return(sim)
})

simdfs22 <- bind_rows(simdf_mod_ts22, .id = "id")
simdfs22$I0 <- param.grid$I0[as.numeric(simdfs22$id)]
simdfs22$R0 <- param.grid$R0[as.numeric(simdfs22$id)]
simdfs22$start_time <- param.grid$start_time[as.numeric(simdfs22$id)]

save(simdfs21, simdfs22, file = "data/sim_start_time.RData")
load("data/sim_start_time.RData")

epi21_day50 <- simdfs_epi_size21 %>%
  filter(start_time==50) %>%
  group_by(R0, I0)

dif21 <- simdfs_epi_size21 %>%
  left_join(epi21_day50, by = c("R0", "I0")) %>%
  mutate(D_dif = D.y - D.x) %>%
  select(I0, R0, D_dif, D.x, D.y, start_time.x)

ggplot(dif21, aes(x = start_time.x, y = I0, z = D_dif)) +
  facet_wrap(~R0) +
  geom_contour_filled()

epi22_day50 <- simdfs_epi_size22 %>%
  filter(start_time==50) %>%
  group_by(R0, I0)

dif22 <- simdfs_epi_size22 %>%
  left_join(epi22_day50, by = c("R0", "I0")) %>%
  mutate(D_dif = D.y - D.x) %>%
  select(I0, R0, D_dif, D.x, D.y, start_time.x)

ggplot(dif22, aes(x = start_time.x, y = I0, z = D_dif)) +
  facet_wrap(~R0) +
  geom_contour_filled()

simdfs_epi_size22 <- simdfs22 %>% filter(time == 120)
ggplot(data = simdfs_epi_size22 %>% filter(R0<5.5), aes(x = start_time, y = I0, z = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_contour_filled(breaks = c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1625, 1650, 1675, 1700, 1703, 1706, 1709, 1712)) +
  labs(x = "Outbreak start time", y = "Initial ratio of infectious individuals", 
       fill = "Final epidemic size\n(cumulative deaths)\nat day 120")
ggsave("figures/2022arriv_start_I0.png")

ggplot(data = simdfs_epi_size22 %>% filter(R0<5.5), aes(x = start_time, y = D)) +
  facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
  geom_line(aes(group = I0, color = I0)) +
  scale_color_viridis_c(name = "Initial ratio of\ninfectious\nindividuals", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/2022_start_D_R0_facets.png", width = 7, height = 5)

ggplot(data = simdfs_epi_size22 %>% filter(R0<5.5 & I0 %in% c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)), aes(x = start_time, y = D)) +
  facet_wrap(~I0, labeller = labeller(I0 = label_both)) +
  geom_line(aes(group = R0, color = R0)) +
  scale_color_viridis_c(name = "R0", option = "plasma") +
  theme_minimal() +
  labs(x = "Outbreak start day (julian date)", y = "Final epidemic size (cumulative deaths) at day 120")
ggsave("figures/2022_start_D_I0_facets.png", width = 7, height = 5)

# Find the start_time that gives maximum D for each R0 and I0 combination
max_death_summary <- simdfs_epi_size22 %>%
  group_by(R0, I0) %>%
  summarise(
    max_death = max(D),
    optimal_start_time = start_time[which.max(D)],
    .groups = 'drop'
  ) %>%
  arrange(R0, I0)

# View the summary
print(max_death_summary)

# Create a heatmap showing optimal start times
ggplot(max_death_summary, aes(x = R0, y = I0, fill = optimal_start_time)) +
  geom_tile() +
    scale_fill_gradient2(name = "Optimal\nStart Time", 
                       low = "blue", mid = "white", high = "red", 
                       midpoint = 50) +
  labs(x = "R0", y = "Initial Infectious Proportion (I0)", 
       title = "Optimal Start Time for Maximum Deaths") +
  theme_minimal()

# Create a heatmap showing maximum deaths achieved
ggplot(max_death_summary, aes(x = R0, y = I0, fill = max_death)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Maximum\nDeaths") +
  labs(x = "R0", y = "Initial Infectious Proportion (I0)", 
       title = "Maximum Deaths Achieved") +
  theme_minimal()

param.grid <- expand.grid(
  start_time = seq(0, 120, by = 1),
  V0 = seq(0, 1, by = 0.1)
)

## simulate vaccination at start time, assuming initial infection = 1
sim_vac22 <- map(1:nrow(param.grid), function(idx) {
  params <- c(beta = 5.39 * 1/12.38, gamma = 1/12.38, mu = 1,
              L = logistic_func[logistic_func$year == 2022, "L"], 
              a = logistic_func[logistic_func$year == 2022, "a"], 
              k = logistic_func[logistic_func$year == 2022, "k"])
  time_vec <- seq(param.grid$start_time[idx], 120)
  init_pop_size <- as.numeric(params["L"] / (1 + exp(-(params["a"] + params["k"] * param.grid$start_time[idx]))))
  init_pop <- c(S = (init_pop_size -1) * (1-param.grid$V0[idx]) , I = 1, R = (init_pop_size - 1) *param.grid$V0[idx], D = 0)
  sim <- as.data.frame(rk4(init_pop, time_vec, mod4, params))
  return(sim)
})
save(sim_vac22, file = "data/sim_vac22.RData")
load("data/sim_vac22.RData")
simdfs22 <- bind_rows(sim_vac22, .id = "id")
simdfs22$V0 <- param.grid$V0[as.numeric(simdfs22$id)]
simdfs22$start_time <- param.grid$start_time[as.numeric(simdfs22$id)]
simdfs22$init_pop <- as.numeric(logistic_func[logistic_func$year == 2022, "L"] / (1 + exp(-(logistic_func[logistic_func$year == 2022, "a"] + logistic_func[logistic_func$year == 2022, "k"] * param.grid$start_time[as.numeric(simdfs22$id)]))))
simdfs22$V0_num <- (simdfs22$init_pop - 1) * param.grid$V0[as.numeric(simdfs22$id)]
simdfs_epi_size22 <- simdfs22 %>% filter(time == 120)

ggplot(data = simdfs_epi_size22, aes(x = start_time, y = V0, z = D))+
  geom_contour_filled()

p0 <- ggplot(data = simdfs_epi_size22 %>% filter(V0 == 1), aes(x = start_time)) +
    # facet_wrap(~V0, label = label_both) +
    geom_line(aes(y = D)) +
    geom_line(aes(y = V0_num), linetype = 2) +
    # geom_line(aes(y = V0_num, group = V0, color = V0), linetype = 2) +
    # geom_line(aes(y = D, group = V0, color = V0)) +
    # scale_color_viridis_c(name = "V0", option = "plasma") +
    scale_y_continuous(name = "Final epidemic size (deaths)", 
                      sec.axis = sec_axis(~., name="Vaccinated number of pelicans")) +
    # labs(x = "Outbreak start time (julian date)") +
    # geom_hline(yintercept = max(simdfs_epi_size22$D)/2, linetype = 3, color = "grey") +
    theme_minimal() +
    theme(axis.title.x = element_blank())

# ggsave("figures/vaccination.png", width = 4, height = 3)

library(gridExtra)
library(cowplot)

p1 <- ggplot(data = simdfs_epi_size22, aes(x = start_time)) +
  geom_line(aes(y = D, group = V0, color = V0)) +
  scale_color_viridis_c(name = "V0", option = "plasma") +
  labs(y = "Final epidemic size (deaths)") +
  theme_minimal() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y.right = element_text())

p2 <- ggplot(data = simdfs_epi_size22, aes(x = start_time)) +
  geom_line(aes(y = V0_num, group = V0, color = V0)) +
  scale_color_viridis_c(name = "Initial vaccinated\nproportion", option = "plasma") +
  labs(x = "Outbreak start time (julian date)", 
       y = "Vaccinated number of birds") +
  theme_minimal() +
  theme(legend.position = "top", axis.title.y.right = element_text())

combined_plot <- plot_grid(p1, p2, ncol = 1, align = "v", 
                          labels = c("A", "B"), label_size = 12)
ggsave("figures/vaccination_V0.png", plot = combined_plot, width = 4, height = 7)

# Calculate the rate of epidemic size decrease per vaccinated individual
vaccination_effectiveness <- simdfs_epi_size22 %>%
  arrange(start_time, V0) %>%
  group_by(start_time) %>%
  summarise(
    # Calculate the slope (rate of change) of D vs V0_num
    rate_D_per_vac = -coef(lm(D ~ V0_num))[2],  # Negative because we want decrease rate
    # Calculate correlation
    correlation = cor(V0_num, D),
    # Calculate effectiveness metrics
    max_deaths = max(D),
    min_deaths = min(D),
    total_reduction = max_deaths - min_deaths,
    max_vac = max(V0_num),
    # Rate per 100 vaccinated
    rate_per_100_vac = rate_D_per_vac * 100,
    .groups = 'drop'
  ) %>%
  filter(!is.na(rate_D_per_vac))
View(vaccination_effectiveness)
# Plot the rate of epidemic reduction vs vaccination coverage
p3 <- ggplot(vaccination_effectiveness, aes(x = start_time, y = rate_per_100_vac)) +
      geom_line() +
      labs(x = "Outbreak start time (julian date)", 
          y = "Deaths prevented\nper 100 vaccinated birds") +
      theme_minimal()


combined_plot2 <- plot_grid(p0, p1, p3, p2, nrow = 2, align = "v", 
                          labels = c("A", "B", "C", "D"), 
                          label_size = 12)
ggsave("figures/vac_comb.png", plot = combined_plot2, width = 8, height = 6)

# fit model
# load("data/2021.RData")
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

param.grid.fit <- expand.grid(
  R0 = seq(1, 15, by = 0.5),
  gamma_inv = seq(1, 15, by = 1),
  I0 = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
  mu = seq(0.5, 1, by = 0.1),
  start_time = c(35, 40, 45)
)

sim_fits <- function(param_grid, data, sir_mod, arriv_func, arriv_year = 2022, arriv_param_confit = NULL){
  simdf_mod = data.frame()
  simdf_mod_timeseries = list()
  simdf_mod_timeseries <- map(1:nrow(param_grid), function(idx) {
    parms <- c(beta = param_grid$R0[idx]/param_grid$gamma_inv[idx], gamma = 1 / param_grid$gamma_inv[idx], mu = param_grid$mu[idx],
               L = arriv_func[arriv_func$year == arriv_year, paste0("L",arriv_param_confit)], 
               a = arriv_func[arriv_func$year == arriv_year, paste0("a", arriv_param_confit)],
               k = arriv_func[arriv_func$year == arriv_year, paste0("k", arriv_param_confit)])
    init_pop_size <- as.numeric(parms["L"] / (1 + exp(-(parms["a"] + parms["k"] * param_grid$start_time[idx]))))
    init_pop <- c(S = init_pop_size * (1 - param_grid$I0[idx]), 
                  I = init_pop_size * param_grid$I0[idx], R = 0, D = 0)
    time_vec <- seq(param_grid$start_time[idx], 120)
    sim <- as.data.frame(rk4(init_pop, time_vec, sir_mod, parms))
    sim_sub <- sim %>% filter(time %in% data$julien_date)
    mse <- sum((sim_sub$D - data$cumdeath)^2) / sum((data$cumdeath)^2)
    sim$mse <- mse
    print(mse)
    return(sim)
  })
  simdf_mod <- bind_rows(simdf_mod_timeseries, .id = "id")
  idx_min <- simdf_mod$id[which.min(simdf_mod$mse)]
  q <- length(param_grid) - 1
  n <- nrow(data)
  simdf_mod_mse_95per <- simdf_mod$mse[which.min(simdf_mod$mse)] + simdf_mod$mse[which.min(simdf_mod$mse)] * q/(n-q)* qf(0.95, n - q, q)
  simdf_mod_conf <- simdf_mod %>% filter(mse < simdf_mod_mse_95per)
  simdf_mod_conf_min <-  simdf_mod_conf %>%
    group_by(time) %>%
    summarise(S = min(S), I = min(I), R = min(R), death = min(D), pop_size = min(pop_size))
  simdf_mod_conf_max <-  simdf_mod_conf %>%
    group_by(time) %>%
    summarise(S = max(S), I = max(I), R = max(R), death = max(D), pop_size = max(pop_size))
  
  p <- ggplot(data = simdf_mod %>% filter(id == idx_min), aes(x = time)) +
    geom_ribbon(aes(ymin = simdf_mod_conf_min$S, ymax = simdf_mod_conf_max$S, fill = "Susceptible"), alpha = 0.2) +
    geom_ribbon(aes(ymin = simdf_mod_conf_min$I, ymax = simdf_mod_conf_max$I, fill = "Infected"), alpha = 0.2) +
    geom_ribbon(aes(ymin = simdf_mod_conf_min$death, ymax = simdf_mod_conf_max$death, fill = "Dead"), alpha = 0.2) +
    geom_ribbon(aes(ymin = simdf_mod_conf_min$pop_size, ymax = simdf_mod_conf_max$pop_size, fill = "Alive"), alpha = 0.2) +
    geom_line(aes(y = S, color = "Susceptible")) +
    geom_line(aes(y = I, color = "Infected")) +
    geom_line(aes(y = D, color = "Dead")) +
    geom_line(aes(y = pop_size, color = "Alive")) +
    labs(x = "Time (days)", y = "Population size", color = "Population", fill = "95% confidence interval") +
    theme_minimal() +
    geom_line(data = data, aes(x = julien_date, y = cumdeath, color = "Observed death"), size = 1) +
    geom_line(data = morta22, aes(x = julien_date, y = TOTAL_DP_ALIVE, color = "Observed alive"), size = 1) +
    scale_color_manual(values = c("black", "darkred", "steelblue", "darkgrey", "pink", "darkgreen")) +
    scale_fill_manual(values = c("black", "darkred", "steelblue", "darkgreen")) +
    annotate("text", x = 90, y = 1000, label = paste("R0:", param_grid[idx_min, "R0"], "\nInfectious period:", 
                                                     round(param_grid[idx_min, "gamma_inv"], 2), "\nMSE:",
                                                     round(simdf_mod$mse[which.min(simdf_mod$mse)], 3), "\nFatality ratio:", 
                                                     param_grid[idx_min, "mu"], "\nI(0):", param_grid[idx_min, "I0"]), size = 5) +
    geom_vline(aes(xintercept = 75), linetype = "dashed", color = "grey") +
    geom_vline(aes(xintercept = 82), linetype = "dashed", color = "grey") +
    geom_vline(aes(xintercept = 100), linetype = "dashed", color = "grey") +
    geom_vline(aes(xintercept = 104), linetype = "dashed", color = "grey")
  return(p)
}

sim_fits(param.grid.fit, morta22, mod4, logistic_func)

## constant infectious proportion of arriving individuals
## all results are NaN
# param.grid2 <- expand.grid(
#   infect_prop = seq(0.1, 1, by = 0.1),
#   R0 = seq(1, 10, by = 0.5),
#   start_time = seq(0, 90, by = 5)
# )
# simdf_mod_ts22_infect_arriv <- map(1:nrow(param.grid2), function(idx) {
#   params <- c(beta = param.grid2$R0[idx] * epi_params[epi_params$year == 2022, "gamma"], 
#               gamma = epi_params[epi_params$year == 2022, "gamma"],
#               mu = epi_params[epi_params$year == 2022, "mu"],
#               infect_prop_in_arriv = param.grid2$infect_prop[idx],
#               L = logistic_func[logistic_func$year == 2022, "L"], 
#               a = logistic_func[logistic_func$year == 2022, "a"], 
#               k = logistic_func[logistic_func$year == 2022, "k"])
#   time_vec <- seq(param.grid2$start_time[idx], 120)
#   init_pop <- c(S = 0, I = 0, R = 0, D = 0)
#   sim <- as.data.frame(rk4(init_pop, time_vec, mod5, params))
#   return(sim)
# })
# 
# simdfs22_ia <- bind_rows(simdf_mod_ts22_infect_arriv, .id = "id")
# simdfs22_ia$infect_prop <- param.grid2$infect_prop[as.numeric(simdfs22_ia$id)]
# simdfs22_ia$R0 <- param.grid2$R0[as.numeric(simdfs22_ia$id)]
# simdfs22_ia$start_time <- param.grid2$start_time[as.numeric(simdfs22_ia$id)]
# 
# simdfs22_ia_epi_size <- simdfs22_ia %>% filter(time == 120)
# ggplot(data = simdfs22_ia_epi_size, aes(x = start_time, y = infect_prop, z = D)) +
#   facet_wrap(~R0, labeller = labeller(R0 = label_both)) +
#   geom_contour_filled() +
#   labs(x = "Outbreak start day (julian date)", y = "Initial infectious individual, I(0)", 
#        fill = "Final epidemic size\n(cumulative deaths)\nat day 120")
