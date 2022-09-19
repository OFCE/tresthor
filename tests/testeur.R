##testeur

library(tresthor)
library(tidyverse)


Sys.setenv("CPATH"="/opt/homebrew/include")
Sys.setenv("LIBRARY_PATH"="/opt/homebrew/lib")
Sys.setenv("PKG_LIBS"="-lsuperlu")

start_solver = as.Date("2010-01-01")
end_solver = as.Date("2020-10-01")

create_model("opale",model_source = "inst/Opale/opale.txt",rcpp= TRUE, rcpp_path = "tests")


coeffs <- readRDS("inst/Opale/coefficients_opale.rds")
data_opale <- readRDS("inst/Opale/donnees_opale.rds")

empty_data <- data_opale %>%
  select(date, all_of(opale@endo_list)) %>%
  mutate_if(is.numeric,~ifelse(date > start_solver, NA,.x)) %>%
  left_join(data_opale %>% select(date,year, all_of(opale@exo_list)))
empty_data<- add_coeffs(coeffs, database = empty_data,pos.coeff.name = 2,pos.coeff.value = 1 )

# Rprof("profil")
solved <- thor_solver(model = opale,first_period = start_solver,last_period = end_solver,
                      main_variable = "td_pib7_ch",main_variable_fx = function(x)(x/lag(x)-1)*100,

                      index_time = "date",rcpp = TRUE,database = empty_data)
# summaryRprof("profil")
data_opale %>% mutate(pib_growth = (td_pib7_ch/lag(td_pib7_ch)-1)*100) %>% filter(date > "2015-01-01") %>%
  select(date,pib_growth)
solved %>% select(year , td_pib7_ch) %>% group_by(year) %>%
  summarise(pib_ann = sum(td_pib7_ch)) %>% filter(year > 2014) %>% mutate(pib_growth = (pib_ann/lag(pib_ann)-1)*100)


t_data[timeref,s_endo]

x_test <- unlist(x_n)
x_test - x_n1
