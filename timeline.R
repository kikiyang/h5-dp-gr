library(ggplot2)
library(lubridate)
load("data/arrival_combined.RData")

demo22 <- data.frame("date" = as.Date(c("2022-01-27", "2022-02-03", "2022-02-07", "2022-02-17", "2022-02-19")), 
                   "event" = c("first arrived pelicans observed", "no pelicans observed",
                               "first continuous increase of arrived pelicans", 
                               "first mortality", "peak arrival rate"))
demo22$year <- 2022

demo21 <- data.frame("date" = as.Date(c("2021-01-05","2021-02-19","2021-02-21")), "event" = c("first nests observed", "first mortality", "peak arrival rate"))
demo21$year <- 2021

genetic22 <- data.frame("date" = as.Date(c("2022-01-20", "2022-02-19", "2022-02-05")), 
                      "event" = c("introduction (95% HPD lower)", 
                                  "introduction (95% HPD higher)", 
                                  "introduction (median)"))
genetic22$year <- 2022

genetic21 <- data.frame("date" = as.Date(c("2020-11-17", "2020-12-28", "2020-12-10")), "event" = c("introduction (95% HPD lower)", 
                                  "introduction (95% HPD higher)", 
                                  "introduction (median)"))
genetic21$year <- 2021

model22 <- data.frame("date" = c(as.Date("2022-02-17")-8.04, 
                                 as.Date("2022-02-17")-16.81, 
                                 as.Date("2022-02-17")-11.62),
                      "event" = c("initial infection (inferred infectious period (95% HPD lower)",
                                  "initial infection (inferred infectious period (95% HPD higher)",
                                  "initial infection (inferred infectious period (median)"))
model22$year <- 2022

model21 <- data.frame("date" = c(as.Date("2021-02-19")-0.11,as.Date("2021-02-19")-0.16, as.Date("2021-02-19")-0.13))
model21$year <- 2021

tl22 <- rbind(demo22, genetic22, model22, demo21, genetic21, model21)

library(ggrepel)
ggplot(tl22, aes(x = date, y = year, label = event)) +
  geom_line() +
  geom_point() +
  geom_text_repel(direction = "y",
                  point.padding = 1.5,
                  hjust = 1.5,
                  box.padding = 1,
                  seed = 55) +
  scale_x_date(name = "", date_breaks = "1 week", date_labels = "%d %B",
               expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_discrete(name = "") +
  theme_minimal()

library(ggplot2)
library(lubridate)
library(ggrepel)

# Separate the data into point events and interval events
demo22 <- data.frame("date" = as.Date(c("2022-01-27", "2022-02-03", "2022-02-07", "2022-02-17", "2022-02-19")), 
                   "event" = c("first arrived pelicans observed", "no pelicans observed",
                               "first continuous increase of arrived pelicans", 
                               "first mortality", "peak arrival rate"))
demo22$year <- 2022

# Create interval data for shaded regions
genetic_interval <- data.frame(
  xmin = as.Date("2022-01-20"),
  xmax = as.Date("2022-02-19"),
  ymin = 2022,
  ymax = 2022.3,
  label = "Introduction period\n(95% HPD)",
  label_x = as.Date("2022-02-05"),
  label_y = 2022
)

model_interval <- data.frame(
  xmin = as.Date("2022-02-17") - 16.81,
  xmax = as.Date("2022-02-17") - 8.04,
  ymin = 2022,
  ymax = 2022.3,
  label = "Initial infection period\n(95% HPD)",
  label_x = as.Date("2022-02-17") - 11.62,
  label_y = 2022
)

# Create the plot
ggplot() +
  # Add shaded boxes for intervals
  geom_rect(data = genetic_interval, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.5, color = NA) +
  geom_rect(data = model_interval, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightcoral", alpha = 0.5, color = NA) +
  
  # Add point events
  geom_line(data = demo22, aes(x = date, y = year)) +
  geom_point(data = demo22, aes(x = date, y = year), size = 2) +
  
  # Add labels for point events
  geom_text_repel(data = demo22, aes(x = date, y = year, label = event),
                  direction = "y",
                  # point.padding = 1.5,
                  # hjust = 1.5,
                  # box.padding = 1,
                  seed = 55) +
  
  # Add median lines for intervals
  geom_vline(xintercept = as.Date("2022-02-05"), linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_vline(xintercept = as.Date("2022-02-17") - 11.62, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_text(aes(x = as.Date("2022-02-05"), y = 2022.32, label = "Virus introduction (phylo)"), 
            color = "blue", angle = 0, vjust = 0, hjust = 1, size = 3) +
  geom_text(aes(x = as.Date("2022-02-17") - 11.62, y = 2022.36, label = "Initial infection (epi)"), 
            color = "red", angle = 0, vjust = 0, hjust = -0.05, size = 3) +
  
  scale_x_date(name = "", date_breaks = "1 week", date_labels = "%d %B",
               expand = expansion(mult = c(0.12, 0.12))) +
  scale_y_continuous(name = "", breaks = 2022, labels = "2022") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(title = "Timeline of HPAI Outbreak Events in 2022")
library(svglite)
ggsave("figures/timeline_2022_intervals.svg", width = 6, height = 3)

load("data/2021.RData")
View(morta21)
View(nest21)
View(arrivRate_comb)

library(ggplot2)
library(lubridate)
library(ggrepel)
library(svglite)

# Point events for 2021
demo21 <- data.frame("date" = as.Date(c("2021-01-05","2021-02-19","2021-02-21")), 
                     "event" = c("first nests observed", "first mortality", "peak arrival rate"))
demo21$year <- 2021

# Create interval data for shaded regions (2021)
genetic_interval_21 <- data.frame(
  xmin = as.Date("2020-11-17"),
  xmax = as.Date("2020-12-28"), 
  ymin = 2021,
  ymax = 2021.1,
  label = "Introduction period\n(95% HPD)",
  label_x = as.Date("2020-12-10"),
  label_y = 2021
)

model_interval_21 <- data.frame(
  xmin = as.Date("2021-02-19") - 0.16,
  xmax = as.Date("2021-02-19") - 0.11,
  ymin = 2021,
  ymax = 2021.1,
  label = "Initial infection period\n(95% HPD)",
  label_x = as.Date("2021-02-19") - 0.13,
  label_y = 2021
)

# Create the 2021 plot
ggplot() +
  # Add shaded boxes for intervals
  geom_rect(data = genetic_interval_21, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.5, color = NA) +
  geom_rect(data = model_interval_21, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightcoral", alpha = 0.5, color = NA) +
  
  # Add point events
  geom_line(data = demo21, aes(x = date, y = year)) +
  geom_point(data = demo21, aes(x = date, y = year), size = 2) +
  
  # Add labels for point events
  geom_text_repel(data = demo21, aes(x = date, y = year, label = event),
                  direction = "y",
                  hjust = 0.5,
                  point.padding = 0.5,
                  box.padding = 0.5,
                  seed = 55) +
  
  # Add median lines for intervals
  geom_vline(xintercept = as.Date("2020-12-10"), linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_vline(xintercept = as.Date("2021-02-19") - 0.13, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_text(aes(x = as.Date("2020-12-10"), y = 2021.1, label = "Virus introduction (phylo)"), 
            color = "blue", angle = 0, vjust = 0, hjust = 0.5, size = 3) +
  geom_text(aes(x = as.Date("2021-02-19") - 0.13, y = 2021.1, label = "Initial infection (epi)"), 
            color = "red", angle = 0, vjust = 0, hjust = 0.5, size = 3) +
  scale_x_date(name = "", date_breaks = "2 weeks", date_labels = "%d %B",
               expand = expansion(mult = c(0.12, 0.12))) +
  theme(axis.ticks.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name = "", breaks = 2021, labels = "2021") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave("figures/timeline_2021_intervals.svg", width = 6, height = 3)
