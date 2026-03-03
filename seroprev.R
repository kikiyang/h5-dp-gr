library(epiR)
# ref: https://cran.r-project.org/web/packages/epiR/vignettes/epiR_diagnostic_tests.html
epi.prev(pos = 2, tested = 20, se = 0.97, sp = 0.98, 
         method = "wilson", units = 1, conf.level = 0.95)

epi.prev(pos = 5, tested = 21, se = 0.97, sp = 0.98, 
         method = "wilson", units = 1, conf.level = 0.95)

epi.prev(pos = 3, tested = 9, se = 0.973, sp = 1,
         method = "wilson", units = 1, conf.level = 0.95)
