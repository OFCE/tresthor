#' Algorithms to determine which equations and variables are part of what
#'
#' @param endogenous_variables vector containing all the declred endogenous variables
#' @param eq_var_matrix matrix of contemporaneous endogenous variables per equations
#' @param decomposition boolean, TRUE if go ahead with decomposition
#'
#' @return list containing all parts infos
#' @keywords internal
#'
#'
decomposing_model<-function(endogenous_variables,eq_var_matrix,decomposition = FALSE){
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

# formula_list<-data.frame(formula=character(),
#                          part=character(),
#                          equation=character(),
#                          stringsAsFactors=FALSE)
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


  #### Epilogue start
  #Keep only equations with the following condition : if a remaining endo appears in only one unique equation, this equation will be in the epilogue.
  ## Remove the equations and test for the others remaining endos
  if(length(remaining_endo)>0){
    occurences<-table(unlist(eqns_mat))

    while(min(occurences)==1){
      epilogue_endo <- c(epilogue_endo, names(occurences[which(occurences==1)]) )

      plop <- as.matrix(eqns_mat)
      plop[!plop %in% epilogue_endo] <- NA
      plop <-as.data.frame(plop)

      epilogue_equations <- c(epilogue_equations , rownames(plop[rowSums(is.na(plop)) != ncol(plop), ]))

      eqns_mat <- eqns_mat[!rownames(eqns_mat) %in% epilogue_equations,]

      occurences<-table(unlist(eqns_mat))
      }
    }
  #### Epilogue ends

  epilogue_endo<-intersect(epilogue_endo,endogenous_variables)
  if(length(epilogue_endo > 0)){epilogue <- TRUE}else{cat("Epilogue block is empty.\n")}

  remaining_endo <- setdiff(remaining_endo,epilogue_endo)
  #### Heart is what's left
  if(length(remaining_endo)>0){
      heart_endo <- remaining_endo
      heart <- TRUE
      heart_equations <- rownames(eqns_mat)
      }else{
      cat("Heart block is empty.\n")}

  ##end of algo
  }

res<- list( prologue=prologue,
            heart=heart,
            epilogue=epilogue,

            prologue_endo=prologue_endo,
            heart_endo=heart_endo,
            epilogue_endo=epilogue_endo,

            prologue_equations=prologue_equations[order(prologue_equations)],
            heart_equations   =heart_equations[order(heart_equations)],
            epilogue_equations=epilogue_equations[order(epilogue_equations)] )

}


