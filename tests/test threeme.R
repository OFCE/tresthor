## ThreeMe Tests

library(tidyverse)
library(tresthor)

classification <-"8x8"

create_model(str_c("threeme",classification),model_source = str_c("tests/threeme_",classification,"_thor.txt"),rcpp = TRUE,rcpp_path = "tests")

data_3me <- readRDS(str_c("tests/data3me_",classification,".rds"))

Rprof("profiler")
data_3me_simulated <- thor_solver(model = get(str_c("threeme",classification)),
                                  first_period =2016,
                                  last_period = 2050,
                                  database = data_3me ,
                                  index_time = "year",
                                  adv_diag = TRUE, rcpp = TRUE)
Rprof(NULL)
summaryRprof("profiler")
saveRDS(data_3me_simulated,str_c("data/threeme",classification,"_simulation_res.rds"))

