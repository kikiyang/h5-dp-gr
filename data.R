rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(nlstools)

# read colony growth data from data_2021_22.xlsx file under data folder
growth <- readxl::read_excel("data/data_2021_22.xlsx", sheet = "2021 subcolony growth")
# pivot growth (colnames are julien dates (days after 2021-01-01), first column is Subcolony)
growth <- growth %>% pivot_longer(cols = -Subcolony, names_to = "julien_date", values_to = "nest_number")
# date of N days after 2021-01-01
growth$julien_date <- as.numeric(gsub("X", "", growth$julien_date)) - 1
growth$date <- as.Date("2021-01-01") + growth$julien_date

# read colony mortality data from data_2021_22.xlsx file under data folder
# AON: Apparently occupied nests
mortality <- readxl::read_excel("data/data_2021_22.xlsx", sheet = "2021 subcolony death")
mortality$Date <- as.Date(mortality$Date, format = "%y-%m-%d")
# Approximating 'many' is 50 and 'a few' is 5
mortality$'Apparently abandoned nests' <- gsub("^many.*", "50", mortality$'Apparently abandoned nests', ignore.case = TRUE)
mortality$'Apparently abandoned nests' <- gsub("^a few*", "5", mortality$'Apparently abandoned nests')
# View(mortality)
mortality <- mortality %>% select(-Area)
# remove non-numeric characters starting from the third column
mortality <- mortality %>% mutate(across(3:ncol(mortality), ~gsub("[^0-9]", "", .)))
mortality <- mortality %>% mutate(across(3:ncol(mortality), as.numeric))
# change NA in mortality data to 0 following the 4 rules below: 1. DEAD ADULTS: No number → 0, 2. DEAD CHICKS: No number → 0, 3. AON Incubated: No number → 0, 4. AON NEW: No number → 0
mortality[is.na(mortality$'Dead adults'), 'Dead adults'] <- 0
mortality[is.na(mortality$'Dead chicks'), 'Dead chicks'] <- 0
mortality[is.na(mortality$'AON incubated, no visible chicks'), 'AON incubated, no visible chicks'] <- 0
mortality[is.na(mortality$'AON new'), 'AON new'] <- 0
# ABANDONED AON: No number → 0 (to the extent it is possible to discriminate abandoned nests-only freshly abandoned nests are visible)
mortality[is.na(mortality$'Apparently abandoned nests'), 'Apparently abandoned nests'] <- 0
# Visible Chicks: No number for those dates prior to the first counted emergencies (chicks) → 0
# Visible Chicks: No number for those dates after emergence of chicks → Not available
# First counted emergencies: 2021-03-26

mortality[is.na(mortality$'Number of visible chicks') & mortality$Date < as.Date("2021-03-26"), 'Number of visible chicks'] <- 0

# Check the consistency between growth and mortality data
# remove (new) and -KORI and /ISL84 from the Subcolony in mortality
mortality$Subcolony <- gsub(" \\(new\\)", "", mortality$Subcolony)
mortality$Subcolony <- gsub("/ISL84", "", mortality$Subcolony)
# uppercase all letters in Subcolony
mortality$Subcolony <- toupper(mortality$Subcolony)


# Add up mortality data for TR46-L, TR46-M+KORI, and TR46-N and rename the Subcolony as TR46
mortality_TR46 <- mortality %>%
    filter(Subcolony == "TR46-L" | Subcolony == "TR46-M+KORI" | Subcolony == "TR46-N") %>%
    group_by(Date) %>%
    summarise(across(c('Dead chicks', 'Dead adults', 'AON incubated, no visible chicks', 'AON new', 'Apparently abandoned nests', 'Number of visible chicks'), sum, na.rm = TRUE)) %>%
    mutate(Subcolony = "TR46")
# Add up mortality data for TR63-KORI and TR63 and rename the Subcolony as TR63
mortality_TR63 <- mortality %>%
    filter(Subcolony == "TR63-KORI" | Subcolony == "TR63") %>%
    group_by(Date) %>%
    summarise(across(c('Dead chicks', 'Dead adults', 'AON incubated, no visible chicks', 'AON new', 'Apparently abandoned nests', 'Number of visible chicks'), sum, na.rm = TRUE)) %>%
    mutate(Subcolony = "TR63")
# Add up mortality data for TR59-KORI and TR59 and rename the Subcolony as TR59
mortality_TR59 <- mortality %>%
    filter(Subcolony == "TR59-KORI" | Subcolony == "TR59") %>%
    group_by(Date) %>%
    summarise(across(c('Dead chicks', 'Dead adults', 'AON incubated, no visible chicks', 'AON new', 'Apparently abandoned nests', 'Number of visible chicks'), sum, na.rm = TRUE)) %>%
    mutate(Subcolony = "TR59")
# Add mortality_TR46 to mortality data and remove TR46-L, TR46-M+KORI, and TR46-N

mortality <- rbind(mortality, mortality_TR46)
mortality <- mortality %>%
    filter(Subcolony != "TR46-L" & Subcolony != "TR46-M+KORI" & Subcolony != "TR46-N" & Subcolony != "TR63-KORI" & Subcolony != "TR63" & Subcolony != "TR59-KORI" & Subcolony != "TR59")
# Add mortality_TR63 to mortality data and remove TR63-KORI and TR63
mortality <- rbind(mortality, mortality_TR63) %>%
    rbind(mortality_TR59)

# cbind growth and mortality data based on Subcolony (if Date is not available in mortality data, add a row)
# Keep all Subcolonies in growth data
colnames(growth) <- c("Subcolony", "julien_date", "nest_number", "Date")

growth_mortality <- merge(growth, mortality, by = c("Subcolony", "Date"), all.x = TRUE, all.y=TRUE) %>%
    select(-julien_date)
growth_mortality <- growth_mortality%>% filter(Subcolony!="TOTAL")

growth %>%
    filter(Subcolony == "TOTAL") %>%
    ggplot(aes(x = julien_date, y = nest_number)) +
    geom_point() +
    geom_line(na.rm = TRUE) +
    labs(x = "Date", y = "Nest number") +
    theme_minimal()

mortality %>%
    group_by(Subcolony, Date) %>%
    summarise(nest_number = cumsum(`AON new`)) %>%
    ungroup() %>%
    group_by(Date) %>%
    summarise(nest_number = sum(nest_number, na.rm = TRUE)) %>%
    ggplot(aes(x = Date, y = nest_number)) +
    geom_point() +
    geom_line(na.rm = TRUE) +
    labs(x = "Date", y = "Nest number") +
    theme_minimal()

## summarise total mortality data for all subcolonies
morta21 <- mortality %>% select(-Subcolony) %>% group_by(Date) %>% 
  summarise(dead_adults = sum(`Dead adults`, na.rm = TRUE), 
            dead_chicks = sum(`Dead chicks`, na.rm = TRUE)) %>%
  mutate(dead_total = dead_adults + dead_chicks, 
         julien_date = as.numeric(Date - as.Date("2021-01-01")),
         cum_death = cumsum(dead_total))

aon_incub <- mortality %>% select(-Subcolony) %>% group_by(Date) %>%
  summarise(aon_incub = sum(`AON incubated, no visible chicks`, na.rm = TRUE),
            aon_new = sum(`AON new`, na.rm = TRUE), abn_nest = sum(`Apparently abandoned nests`, na.rm = TRUE)) %>%
  mutate(julien_date = as.numeric(Date - as.Date("2021-01-01")))


## fit a logistic growth curve for nest number (*2)
nest21 <- growth %>% filter(Subcolony=="TOTAL")
nest21 <- nest21 %>% select(-Subcolony)
## add a new column: nest_number_adjusted, when julien_date > 49, nest_number_adjusted = aon_incub; otherwise, nest_number_adjusted = nest_number
## adjusting vantage point ob (*1.5) to drone ob
nest21 <- left_join(nest21, aon_incub, by = "julien_date")
nest21 <- nest21 %>% mutate(nest_number_adjusted = ifelse(is.na(Date.y), nest_number * 1.5, aon_incub + abn_nest))
nest21$est_alive <- nest21$nest_number_adjusted * 2
nest21 <- nest21 %>% filter(julien_date < 120)

save(nest21, morta21, file = "data/2021.RData")
load("data/2021.RData")

## fitting a logistic growth model
logist_pars21 <- coef(lm(qlogis(nest21$est_alive/2500) ~ julien_date, data = nest21))
logist_func21 <- nls(est_alive~L/(1+exp(-(a+k*julien_date))), start = list(L=2500,a=logist_pars21[1],k=logist_pars21[2]), 
                   data = nest21, trace = TRUE)
summary(logist_func21)
L21 <- coef(logist_func21)[1]
a21 <- coef(logist_func21)[2]
k21 <- coef(logist_func21)[3]
x21 <- c(0:max(nest21$julien_date))
y21 <- L21/(1+exp(-(a21+k21*x21)))

logist_bs21 <- nlsBoot(logist_func21)
ci21 <- logist_bs21$bootCI
L_low21 <- ci21["L", "2.5%"]
L_high21 <- ci21["L", "97.5%"]
a_low21 <- ci21["a", "2.5%"]
a_high21 <- ci21["a", "97.5%"]
k_low21 <- ci21["k", "2.5%"]
k_high21 <- ci21["k", "97.5%"]
y_low21 <- L_low21 / (1 + exp(-(a_low21 + k_low21 * x21)))
y_high21 <- L_high21 / (1 + exp(-(a_high21 + k_high21 * x21)))
predict21 <- data.frame(x21, y21, y_low21, y_high21)

f <- expression(L / (1 + exp(-(a + k * x))))
fd <- D(f, "x")
arriv21 <- data.frame(julien_date = x21, arriv_rate = eval(fd, list(L = L21, a = a21, k = k21, x = x21)))

## 2020 data
morta20 <- readxl::read_excel("data/data_2021_22.xlsx", sheet = "2020 drone")
morta20$Date <- as.Date(morta20$Date, format = "%y-%m-%d")
morta20 <- morta20 %>% mutate(across(4:ncol(morta20), ~gsub("~", "", .)))
morta20 <- morta20 %>% mutate(across(4:ncol(morta20), ~gsub("UNKNOWN", 9999, .)))
morta20 <- morta20 %>% mutate(across(4:ncol(morta20), as.numeric))
morta20[is.na(morta20$'Dead adults'), 'Dead adults'] <- 0
morta20[is.na(morta20$'Dead chicks'), 'Dead chicks'] <- 0
morta20[is.na(morta20$'AON incubated, no visible chicks'), 'AON incubated, no visible chicks'] <- 0
morta20[is.na(morta20$'AON new'), 'AON new'] <- 0
morta20[is.na(morta20$'Apparently abandoned nests'), 'Apparently abandoned nests'] <- 0
morta20[is.na(morta20$'Number of visible chicks') & morta20$Date < as.Date("2021-04-10"), 'Number of visible chicks'] <- 0
morta20[morta20 == 9999] <- NA

# calculate present nest number (sum of AON new before the Date - Apparently abandoned nests)
morta20 <- morta20 %>%
            group_by(Subcolony) %>%
            mutate(nest_number = cumsum(`AON new`) - cumsum(`Apparently abandoned nests`)) %>%
            ungroup()

morta20_sumnest <- morta20 %>%
    group_by(Date) %>%
    summarise(nest_number = sum(nest_number, na.rm = TRUE))

growth20_vintage <- readxl::read_excel("data/data_2021_22.xlsx", sheet = "2020 nest vintage")
growth20 <- growth20_vintage %>% pivot_longer(cols = -Date, names_to = "subcolony", values_to = "nest_number")

growth20$Date <- as.Date(growth20$Date, format = "%y-%m-%d")
growth20_sum <- growth20 %>%
    group_by(Date) %>%
    filter(Date < as.Date("2020-03-06")) %>%
    summarise(nest_number = sum(nest_number, na.rm = TRUE)*701/523*1.6)
nest20 <- rbind(morta20_sumnest,growth20_sum)

# julien date of Date
nest20$julien_date <- as.numeric(nest20$Date - as.Date("2020-01-01"))

nest20 <- nest20[nest20$Date < as.Date("2020-07-01"), ]
nest20$ind_no <- 2 * nest20$nest_number
save(nest20,file = "data/2020.RData")
load("data/2020.RData")

logist_pars20 <- coef(lm(qlogis(ind_no/2800) ~ julien_date, data = nest20))
logist_func20 <- nls(ind_no~L/(1+exp(-(a+k*julien_date))), start = list(L=2800,a=logist_pars20[1],k=logist_pars20[2]), 
                data = nest20, trace = TRUE)
summary(logist_func20)
L20 <- coef(logist_func20)[1]
a20 <- coef(logist_func20)[2]
k20 <- coef(logist_func20)[3]
x20 <- c(0:max(nest20$julien_date))
y20 <- L20/(1+exp(-(a20+k20*x20)))

logist_func_bs20 <- nlsBoot(logist_func20)
ci20 <- logist_func_bs20$bootCI
L_low20 <- ci20["L", "2.5%"]
L_high20 <- ci20["L", "97.5%"]
a_low20 <- ci20["a", "2.5%"]
a_high20 <- ci20["a", "97.5%"]
k_low20 <- ci20["k", "2.5%"]
k_high20 <- ci20["k", "97.5%"]
y_low20 <- L_low20/(1+exp(-(a_low20+k_low20*x20)))
y_high20 <- L_high20/(1+exp(-(a_high20+k_high20*x20)))
predict20 <- data.frame(x20, y20, y_low20, y_high20)

arriv20 <- data.frame(julien_date = x20, arriv_rate = eval(fd, list(L = L20, a = a20, k = k20, x = x20)))

# read 2022 data
morta22 <- readxl::read_excel("data/data_2021_22.xlsx", sheet = "2022 subcolony death")
morta22$DATE <- as.Date(morta22$DATE, format = "%y-%m-%d")
morta22 <- morta22 %>% mutate(nest_number = cumsum(TOTAL_NESTS))
morta22$julien_date <- as.numeric(morta22$DATE - as.Date("2022-01-01"))
# NA in morta22$TOTAL_DEAD_DP to 0
morta22[is.na(morta22$TOTAL_DEAD_DP), 'TOTAL_DEAD_DP'] <- 0
# adjust vantage point data of total pelicans alive by multiplying by 1.49 and rounding to the nearest integer
morta22[morta22$Observation_method == "Vintage", 'TOTAL_DP_ALIVE'] <- 
    round(morta22[morta22$Observation_method == "Vintage", 'TOTAL_DP_ALIVE'] 
          * (1089 +259)/882) # From comparison of data on Feb 25 between drone images and vantage point
morta22 <- morta22 %>%
    mutate(cumdeath = cumsum(TOTAL_DEAD_DP)) %>%
    ungroup() %>%
    mutate(theor_pop = cumdeath + TOTAL_DP_ALIVE)

morta22[morta22$Observation_method=="NO","TOTAL_DEAD_DP"] <- NA
morta22[morta22$DATE=="2022-04-11","TOTAL_DEAD_DP"] <- NA

save(morta22, file = "data/2022.RData")
load("data/2022.RData")

morta22_simp <- morta22 %>% filter(!is.na(TOTAL_DEAD_DP)) %>% select(DATE, TOTAL_DEAD_DP, julien_date, cumdeath) %>%
                mutate(death_adults = NA, death_chicks = NA)
colnames(morta22_simp) <- c("Date", "dead_total", "julien_date", "cum_death", "dead_adults", "dead_chicks")
morta21_22 <- rbind(morta22_simp, morta21) %>%
    mutate(year = lubridate::year(Date))

p21 <- ggplot(morta21, aes(x=julien_date)) +
  geom_point(aes(y = cum_death/2.5, shape = "Cumulative deaths"), size = 3) +
  geom_line(aes(y = cum_death/2.5, linetype = "Cumulative deaths")) +
  geom_point(aes(y = dead_chicks, shape = "Chicks (New death)"), size = 3) +
  geom_line(aes(y = dead_chicks, linetype = "Chicks (New death)")) +
  geom_point(aes(y = dead_adults, shape = "Adults (New death)"), size = 3) +
  geom_line(aes(y = dead_adults, linetype = "Adults (New death)")) +
  scale_y_continuous(sec.axis = sec_axis(~.*2.5, name = "Cumulative deaths"), limits = c(0, 700)) +
  scale_shape_manual(values = c("Cumulative deaths" = 16, "Chicks (New death)" = 2, "Adults (New death)" = 3)) +
  scale_linetype_manual(values = c("Cumulative deaths" = "solid", "Chicks (New death)" = "dashed", "Adults (New death)" = "dotted")) +
  labs(x = "Julian date", y = "New deaths", shape = "Category", linetype = "Category") +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.6),  # x=0.8 (right), y=0.5 (middle)
        legend.background = element_rect(fill = "transparent", color = "white", size = 0.5),
        legend.margin = margin(6, 6, 6, 6), plot.title = element_text(hjust = 0.5)) +
  guides(shape = guide_legend(title = NULL), 
         linetype = guide_legend(title = NULL)) +
  xlim(0, 210) +
  ggtitle("2021 outbreak (H5N8, Clade 2.3.4.4b, EA-2020-A)")

p22 <- ggplot(morta22_simp, aes(x=julien_date)) +
  geom_point(aes(y = cum_death/2.5, shape = "Cumulative deaths"), size = 3) +
  geom_line(aes(y = cum_death/2.5, linetype = "Cumulative deaths")) +
  geom_point(aes(y = dead_total, shape = "New deaths"), size = 3) +
  geom_line(aes(y = dead_total, linetype = "New deaths")) +
  scale_y_continuous(sec.axis = sec_axis(~.*2.5, name = "Cumulative deaths"), limits = c(0, 700)) +
  scale_shape_manual(values = c("Cumulative deaths" = 16, "New deaths" = 1)) +
  scale_linetype_manual(values = c("Cumulative deaths" = "solid", "New deaths" = "dashed")) +
  labs(x = "Julian date", y = "New deaths", shape = "Category", linetype = "Category") +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.6), 
        legend.background = element_rect(fill = "transparent", color = "transparent", size = 0.5),
        legend.margin = margin(6, 6, 6, 6), plot.title = element_text(hjust = 0.5)) +
  guides(shape = guide_legend(title = NULL), 
         linetype = guide_legend(title = NULL))+
  xlim(0, 210) +
  ggtitle("2022 outbreak (H5N1, Clade 2.3.4.4b, EA-2021-AB)")

p21+p22
ggsave("figures/2021_2022_death_scaled.png", height = 6, width = 10)

load("data/2022.RData")
arriv22 <- morta22 %>% filter(!is.na(theor_pop) & julien_date >37)
# arriv22 <- morta22 %>% filter(Corpse_removal == "NO") %>% filter(!is.na(theor_pop)) %>% filter(theor_pop!=0)
logist_pars22 <- coef(lm(qlogis(theor_pop/2500) ~ julien_date, data = arriv22))
logist_func22 <- nls(theor_pop~L/(1+exp(-(a+k*julien_date))), start = list(L=2500,a=logist_pars22[1],k=logist_pars22[2]), 
                data = arriv22, trace = TRUE)
summary(logist_func22)
L22 <- coef(logist_func22)[1]
a22 <- coef(logist_func22)[2]
k22 <- coef(logist_func22)[3]
x22 <- c(0:max(arriv22$julien_date))
y22 <- L22/(1+exp(-(a22+k22*x22)))

# bootstrap the data and generate confidence intervals
logist_func_bs22 <- nlsBoot(logist_func22)
L_low22 <- logist_func_bs22$bootCI["L", "2.5%"]
L_high22 <- logist_func_bs22$bootCI["L", "97.5%"]
a_low22 <- logist_func_bs22$bootCI["a", "2.5%"]
a_high22 <- logist_func_bs22$bootCI["a", "97.5%"]
k_low22 <- logist_func_bs22$bootCI["k", "2.5%"]
k_high22 <- logist_func_bs22$bootCI["k", "97.5%"]
y_low22 <- L_low22/(1+exp(-(a_low22+k_low22*x22)))
y_high22 <- L_high22/(1+exp(-(a_high22+k_high22*x22)))
predict22 <- data.frame(x22, y22, y_low22, y_high22)

growth_rate <- data.frame(julien_date = x22, growth_rate = eval(fd, list(L = L22, a = a22, k = k22, x = x22)))

# # 2020, 2021, 2022 arrival rate combined
arriv22 <- arriv22 %>% select(DATE, nest_number, julien_date, theor_pop)
colnames(arriv22) <- c("Date", "nest_number", "julien_date", "est_alive")
colnames(nest20) <- c("Date", "nest_number", "julien_date", "est_alive")
nest21 <- nest21 %>% select(-Date.y, -nest_number, -aon_new, -aon_incub, -abn_nest)
colnames(nest21) <- c("julien_date", "Date", "nest_number", "est_alive")
nest21 <- nest21 %>% arrange(Date, nest_number, julien_date, est_alive)
arriv_comb <- rbind(arriv22, nest21, nest20)
arriv_comb$year <- lubridate::year(arriv_comb$Date)

predict20$year <- "2020"
predict21$year <- "2021"
predict22$year <- "2022"
colnames(predict20) <- c("x", "y", "y_low", "y_high", "year")
colnames(predict21) <- c("x", "y", "y_low", "y_high", "year")
colnames(predict22) <- c("x", "y", "y_low", "y_high", "year")
predict_comb <- rbind(predict20, predict21, predict22)

growth_rate$year <- "2022"
arriv20$year <- "2020"
arriv21$year <- "2021"
colnames(growth_rate) <- c("julien_date", "arriv_rate", "year")
arrivRate_comb <- rbind(arriv20, arriv21, growth_rate)

logistic_func <- data.frame(
    year = c("2020", "2021", "2022"),
    L = c(L20, L21, L22),
    a = c(a20, a21, a22),
    k = c(k20, k21, k22)
)

save(arriv_comb, predict_comb, arrivRate_comb, logistic_func, file = "data/arrival_combined.RData")
load("data/arrival_combined.RData")
ggplot(arriv_comb)+
    facet_wrap(~year, ncol = 3)+
    geom_point(aes(x=julien_date, y=est_alive), color='blue',size=5)+
    labs(x='Julian date',y='Arrivals')+
    geom_line(data = predict_comb, aes(x=x, y=y), size=1) +
    geom_line(data = predict_comb, aes(x=x, y=y_low), linetype = "dashed", color="lightblue") +
    geom_line(data = predict_comb, aes(x=x, y=y_high), linetype = "dashed", color="lightblue") +    
    geom_ribbon(data = predict_comb, aes(x=x, ymin=y_low, ymax=y_high), fill = "lightblue", alpha = 0.25) +
    geom_line(data=arrivRate_comb, aes(x = julien_date, y = 20*arriv_rate), color = "red", size = 1) +
    scale_y_continuous(sec.axis = sec_axis(~./20, name = "Arrival rate (ind./day)")) +
    theme_bw() +
    theme(axis.text=element_text(size=18),axis.title=element_text(size=24), strip.text = element_text(size = 18))
ggsave("figures/arrival_allYears.png", width = 15, height = 6)
