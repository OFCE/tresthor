library(tresthor)
library(tidyverse)

## Opale

start_solver = as.Date("1983-01-01")
end_solver = as.Date("2020-10-01")

create_model("opale",model_source =system.file("Opale/opale.txt",package = "tresthor"),rcpp = TRUE,rcpp_path = "tests")
# save_model(opale,folder_path = "models/")
load_model(file = "models/opale.rds")

coeffs <- readRDS(system.file("Opale/coefficients_opale.rds",package = "tresthor"))
data_opale <- readRDS(system.file("Opale/donnees_opale.rds",package = "tresthor"))
empty_data <- data_opale %>%
  select(date, all_of(opale@endo_list)) %>%
  mutate_if(is.numeric,~ifelse(date > "1983-01-01", NA,.x)) %>%
  left_join(data_opale %>% select(date,year, all_of(opale@exo_list)))
empty_data<- add_coeffs(coeffs, database = empty_data,pos.coeff.name = 2,pos.coeff.value = 1 )

Rprof("profil")
solved <- thor_solver(model = opale,first_period = start_solver,last_period = end_solver,
                      main_variable = "td_pib7_ch",main_variable_fx = function(x)(x/lag(x)-1)*100,index_time = "date",rcpp = TRUE,database = empty_data)
Rprof(NULL)
summaryRprof("profil")


## Mesange
start_solver = "1995Q1"
end_solver = "2100Q4"
## Voir .dynlib
Rprof("profil")
create_model("mesange",model_source = "tests/mes17.txt",rcpp = TRUE,rcpp_path = "tests")
create_model("mesange_superlu",model_source = "tests/mes17.txt",rcpp = TRUE,rcpp_path = "tests")
Rprof(NULL)
summaryRprof("profil")
save_model(mesange,folder_path = "tests/")
load_model(file = "tests/mesange.rds")

data_mesange <- readRDS("data/mesange_data.rds")
coeff_mes <- read.csv("data/coeff_Mesange_2017.csv",sep = ";")

empty_data <- data_mesange %>%
  select(date, all_of(mesange@endo_list)) %>%
  mutate_if(is.numeric,~ifelse(date > "1995Q1", NA,.x)) %>%
  left_join(data_mesange %>% select(date,year, all_of(mesange@exo_list)))


empty_data<- add_coeffs(coeff_mes, database = empty_data )

Rprof("profil")
solved <- thor_solver(model = mesange,
                      first_period = start_solver,
                      last_period = end_solver,
                      index_time = "date",
                      rcpp = TRUE,
                      database = empty_data)
Rprof(NULL)
summaryRprof("profil")


