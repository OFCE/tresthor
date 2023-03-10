### Faster decomposition algo
endogenous_variables <- endo
eq_var_matrix <-eqns
decomposition = TRUE
objective <-readRDS("tests/threemedecomposed")

eqns_as_list<-purrr::map(as.data.frame(t(eqns)),~unique(stats::na.omit(.x)))

eqns_list <- eqns_as_list %>% purrr::map(~list(endo = .x, n_endo = length(.x)))

prologue_left <-

### 0.A initialising the outputs

prologue<-FALSE
heart<-FALSE
epilogue<-FALSE

prologue_endo <- as.character(vector())
heart_endo <- as.character(vector())
epilogue_endo <- as.character(vector())

prologue_equations <- as.character(vector())
heart_equations <- as.character(vector())
epilogue_equations <- as.character(vector())

### 0.B setup
remaining_endo <- endogenous_variables
equations_vector <- rownames(eq_var_matrix)
eqns_mat <- as.data.frame(eq_var_matrix)


### 1 No algorithm used : only one block is generated
if(decomposition == FALSE){

  heart <- TRUE
  heart_endo <- remaining_endo
  heart_equations <- equations_vector

}else{

  ### 2 Algorithm to decompose in three blocks

  #### Prologue start
  #### Prologue : All equations that contain only one endogenous variable OR only endogenous already determined to be in the prologue
  eqns_mat$nbendo <- 0

  ##count how many endogenous are left in the equations (-1 cause one column is the nb endo counter)
  eqns_mat$nbendo <- apply(eqns_mat, 1 , function(y) length(which(is.na(y)==FALSE)) )-1


  while(min(eqns_mat$nbendo) == 1){

    unique_endo <- eqns_mat %>% dplyr::filter(nbendo == 1) %>% dplyr::select(-nbendo)

    prologue_equations <- c(prologue_equations,rownames(unique_endo))
    prologue_endo <- c(prologue_endo, stats::na.omit(unique(as.vector(as.matrix(unique_endo))) ))

    eqns_mat <- as.matrix(eqns_mat)
    eqns_mat[eqns_mat %in% prologue_endo] <- NA
    eqns_mat <- as.data.frame(eqns_mat)

    eqns_mat$nbendo <- apply(eqns_mat, 1 , function(y) length(which(is.na(y)==FALSE)) )-1
    eqns_mat <- eqns_mat %>% dplyr::filter(nbendo > 0)
    if(length(eqns_mat$nbendo)==0){break}
  }

  if(length(prologue_endo > 0)){prologue <- TRUE}else{cat("Prologue block is empty.\n")}

  eqns_mat$nbendo <- NULL

  #### Prologue end

  remaining_endo <- setdiff(remaining_endo,prologue_endo)


}
