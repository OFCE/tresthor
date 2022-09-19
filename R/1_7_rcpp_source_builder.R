##Creating the jacobian and equations functions for rcpp
# path_jacobian_R= "temp_paprfn/p_j_f.R"
# var_order = all_model_variables
# index_jacobienne =prologue_jacobian
# cpp_file = NA
# path_f_x = "temp_paprfn/p_e_f.R"

#' Transforming the jacobian R source code into Rcpparmadillo code
#'
#' @param path_jacobian_R where and which jacobian
#' @param var_order all variables in alphabetical order
#' @param index_jacobienne symbolic jacobian
#' @param cpp_file name and location of the file.
#' @author Charlotte Nudelmann (DG Tresor)
#' @importFrom purrr map_chr
#' @import stringr
#' @return cpp script
#' @keywords internal
#'
Jacobian_Rcpp <- function(path_jacobian_R,
                          var_order = var_list_sorted,
                          index_jacobienne,
                          cpp_file = NA){

  #Copie du fichier R en txt puis chargement du fichier txt
  file.copy(path_jacobian_R, "jacobian.txt")
  jacobian <- scan("jacobian.txt",character(),sep="\n",quiet = TRUE)
  file.remove("jacobian.txt")

  firstFour <- jacobian[1:4]

  last <- jacobian[length(jacobian)]

  #traiter les premieres lignes
  firstTreated <- character()
  firstTreated[1] <- paste("mat ", gsub("<.+$","",firstFour[1]),"(arma::mat & M, int & t){",sep="")
  firstTreated[2] <- paste("mat Jacobian_n(",
                           substring(stringr::str_extract(string = firstFour[2], pattern = "= \\d+")[[1]],3),
                           ",",
                           substring(stringr::str_extract(string = firstFour[2], pattern = "= \\d+")[[1]],3),
                           ",arma::fill::zeros);", sep="")
  positionsAvecUn <- which(index_jacobienne=='1')
  firstTreated[3] <- paste("Jacobian_n.elem( arma::uvec(\"",
                           paste(positionsAvecUn-1,collapse = ' '),
                           "\") ) = arma::vec(",
                           length(positionsAvecUn),
                           ",arma::fill::ones);",
                           sep = "")
  positionsAvecMoinsUn <- which(index_jacobienne=='-1')
  firstTreated[4] <- paste("Jacobian_n.elem( arma::uvec(\"",
                           paste(positionsAvecMoinsUn-1,collapse = ' '),
                           "\") ) = arma::vec(",
                           length(positionsAvecMoinsUn),
                           ",arma::fill::ones) - 2;",
                           sep = "")

  #traiter la derniere ligne
  lastTreated <- c("return(Jacobian_n);","}")

  #
  jacobian <- jacobian[5:(length(jacobian)-1)]

  #suppression des espaces
  jacobian <- as.character(vapply(jacobian,function(x)gsub(" ","",x), character(1)))

  #traiter les elements e1
  index_a_traiter <- grep("<-\\(\\{", jacobian)
  #Amelioration, si jamais il y a du e2, e3, e4 etc.
  #e_index <- 1
  # while (any(grepl(paste("e",e_index,sep=""), jacobian))) {
  #
  #
  #   e_index = e_index + 1
  # }
  for(i in index_a_traiter){
    e1 <- gsub("e1<-","",jacobian[i+1])
    jacobian[i+2] <-gsub("e1",paste("(",e1,")", sep =""),jacobian[i+2])
    jacobian[i] <- gsub("\\{|\\}","",paste(jacobian[i],jacobian[i+2],jacobian[i+3],sep=""))
    jacobian[i+1] <- NA
    jacobian[i+2] <- NA
    jacobian[i+3] <- NA
  }
  jacobian <- jacobian[!is.na(jacobian)]
  toto<-as.data.frame(jacobian)
  toto <- toto[!is.na(toto)]
  #integrer la date t dans t_data[,variable] et suppression du [t] en fin de ligne
  jacobian <- as.character(vapply(jacobian,function(x)gsub("t\\_data\\[,","t\\_data\\[t,",x), character(1)))
  jacobian <- as.character(vapply(jacobian,function(x)gsub("\\[t\\]","",x), character(1)))
  for (i in 1:length(jacobian)){
    jacobian[i]<-paste(jacobian[i],";",sep="")
  }
  #jacobian <- as.character(vapply(jacobian,function(x)gsub("\\[t\\]",";",x), character(1)))

  #on supprime la fonction lag
  jacobian <- as.character(vapply(jacobian,function(x)gsub("lag\\(t_data\\[t(,'[[:alnum:]|_]+')\\],(\\d)\\)","t_data\\[t-\\2\\1\\]",x), character(1)))
  jacobian <- as.character(vapply(jacobian,function(x)gsub("lag\\(t_data\\[t(,[[:alnum:]|_]+)\\],(t_data\\[t(,[[:alnum:]|_]+)\\]\\+?\\d?)\\)","t_data\\[t-(\\2)\\1\\]",x), character(1)))




  #remplacer la variable par son index
  substituteVariable <- function(text, var_order){
    index_vars<-purrr::set_names(order(var_order)-1,var_order)
    varInTextReplace <- stringr::str_extract_all(string = text, pattern = "'([:alnum:]|_)+'")[[1]]
    for (i in varInTextReplace){
      clean_i <- gsub("'","",i)
      text <- gsub(i,index_vars[clean_i], text) #-2 car c++ commence a 0 et R a 1 et qu'on enleve la premiere colonne (date) dans C++
    }
    return <- text
  }

  jacobian<- purrr::map_chr(jacobian,~substituteVariable(.x,var_order))
  #remplacer dataframe[] par submat de RccpArmadillo
  jacobian <- as.character(vapply(jacobian,function(x)gsub("t_data\\[(t-?\\d?-?\\d?)(,\\d+)\\]","M\\(\\1\\2\\)",x), character(1)))
  jacobian <- as.character(vapply(jacobian,function(x)gsub("t_data\\[(t-?\\d?-?\\d?-\\(M\\(t-?\\d?,\\d+\\)\\))(,\\d+)\\]","M\\(\\1\\2\\)",x), character(1)))

  jacobian <- as.character(vapply(jacobian,function(x)gsub("Jacobian_n\\[(\\d+)(,\\d+)\\]","Jacobian_n\\(\\1-1\\2-1\\)",x), character(1)))


  #remplacer <-  par =
  jacobian <- as.character(vapply(jacobian,function(x)gsub("<-","=",x), character(1)))

  jacobian_f <- c(firstTreated, jacobian, lastTreated, "")


  #remplacer les ^n par pow(...,n)
  for (i_text in 1:length(jacobian_f)){
    while(grepl("\\^\\d",jacobian_f[i_text])){ #s'il y a une puissance
      pos_parenthese_ouverte <- integer()
      pos_power <- gregexpr("\\^\\d",jacobian_f[i_text])[[1]][1]-2# on se positionne avant de fermer la parenthese de power
      #tant qu'il reste du newdiff a traiter
      #Pour chaque newdiff on veut reperer la parenthese fermante
      nb_parenthese_ouverte <- -1
      index_char <- pos_power
      while (nb_parenthese_ouverte < 0){
        nb_parenthese_ouverte <- nb_parenthese_ouverte +
          (substring(jacobian_f[i_text],index_char,index_char)=="(") -
          (substring(jacobian_f[i_text],index_char,index_char)==")")
        index_char <- index_char - 1
      }
      pos_parenthese_fermee <- index_char+1
      if(substring(jacobian_f[i_text],pos_parenthese_fermee-1,pos_parenthese_fermee-1)=="M"){
        pos_parenthese_fermee=pos_parenthese_fermee-1
      }
      power <- substring(jacobian_f[i_text],pos_power+3,pos_power+3)
      jacobian_f[i_text]<-gsub(gsub("(\\[)","\\\\\\1",gsub("(\\])","\\\\\\1",gsub("(\\))","\\\\\\1",gsub("(\\()","\\\\\\1",gsub("(\\+)","\\\\\\1",gsub("(\\*)","\\\\\\1",gsub("(\\^\\d+)","\\\\\\1",
                                                                                                                                                                              substring(jacobian_f[i_text],pos_parenthese_fermee,pos_power+3)))))))),
                               paste("pow(",
                                     substring(jacobian_f[i_text],pos_parenthese_fermee,pos_power+1),
                                     ",",power,")",
                                     sep=""),
                               jacobian_f[i_text])
    }
  }

  if(!is.na(cpp_file)){
    if(file.exists(paste(gsub("<.+$","",firstFour[1]),".cpp",sep=""))){file.remove(paste(gsub("<.+$","",firstFour[1]),".cpp",sep=""))}
    cat(jacobian_f,file=cpp_file,sep="\n",append = TRUE)
  }else{
    if(file.exists(paste(gsub("<.+$","",firstFour[1]),".cpp",sep=""))){file.remove(paste(gsub("<.+$","",firstFour[1]),".cpp",sep=""))}
    cat(jacobian_f,file=paste(gsub("<.+$","",firstFour[1]),".cpp",sep=""),sep="\n",append = FALSE)
  }
}

#' Transforming the equations R source code into Rcpparmadillo code
#'
#' @param path_f_x where and which equation source file
#' @param var_order all model variables in alphabetical order
#' @param cpp_file name of cpp file where to write
#' @importFrom purrr map
#' @import stringr
#' @author Charlotte Nudelmann (DG Tresor)
#' @keywords internal
#' @return cpp script
f_x_Rcpp <- function(path_f_x, var_order, cpp_file = NA){

  file.copy(path_f_x, "f_x.txt")
  f_x <- unique(scan("f_x.txt",character(),sep="\n",quiet=TRUE))
  file.remove("f_x.txt")

  #traiter les deux premieres lignes
  first <- f_x[1:2]

  firstTreated <- character()
  firstTreated[1] <- paste("arma::vec ", gsub("<.+$","",first[1]),"(arma::mat& M, int& t){",sep="")
  firstTreated[2] <- paste("arma::vec f_x_n(",
                           substring(stringr::str_extract(string = first[2], pattern = ",\\d+\\)")[[1]],2,nchar(stringr::str_extract(string = first[2], pattern = ",\\d+\\)")[[1]])-1),
                           ");", sep="")


  #traiter la derniere ligne
  lastTreated <- c("return(f_x_n);","}")

  #Traitement du corps
  f_x <- f_x[3:(length(f_x)-1)]

  #suppression des espaces
  f_x <- as.character(vapply(f_x,function(x)gsub(" ","",x), character(1)))

  #remplacer la variable par son index
  substituteVariable <- function(text, varList){
    varInTextReplace <- stringr::str_extract_all(string = text, pattern = "'([:alnum:]|_)+'")[[1]]
    for (i in varInTextReplace){
      text <- gsub(i,grep(paste("^",substring(i,2,nchar(i)-1),"$",sep=""),varList)-1, text) #-2 car c++ commence a 0 et R a 1 et qu'on enleve la premiere colonne (date) dans C++
    }
    return <- text
  }
  f_x<- purrr::map_chr(f_x,~substituteVariable(.x,var_order))


  #integrer la date t dans t_data[,variable] et suppression du [t] en fin de ligne
  f_x <- as.character(vapply(f_x,function(x)gsub("t\\_data\\[,","t\\_data\\[t,",x), character(1)))
  f_x <- as.character(vapply(f_x,function(x)gsub("\\[t\\]",";",x), character(1)))

  #on supprime la fonction lag
  f_x <- as.character(vapply(f_x,function(x)gsub("lag\\(t_data\\[t(,[[:alnum:]|_]+)\\],(\\d)\\)","t_data\\[t-\\2\\1\\]",x), character(1)))
  f_x <- as.character(vapply(f_x,function(x)gsub("lag\\(t_data\\[t(,[[:alnum:]|_]+)\\],(t_data\\[t(,[[:alnum:]|_]+)\\]\\+?\\d?)\\)","t_data\\[t-(\\2)\\1\\]",x), character(1)))


  #on supprime les fonctions new_diff
  #1 reperer les newdiff
  for (i_text in 1:length(f_x)){
    while(grepl("newdiff",f_x[i_text])){ #s'il y a un newdiff
      pos_parenthese_fermee <- integer()
      pos_newdiff <- gregexpr("newdiff\\(",f_x[i_text])[[1]][1]+8# on se positionne au debut de l'argument de newdiff
      #tant qu'il reste du newdiff a traiter
      #Pour chaque newdiff on veut reperer la parenthese fermante
      nb_parenthese_ouverte <- 1
      index_char <- pos_newdiff
      while (nb_parenthese_ouverte > 0){
        nb_parenthese_ouverte <- nb_parenthese_ouverte +
          (substring(f_x[i_text],index_char,index_char)=="(") -
          (substring(f_x[i_text],index_char,index_char)==")")
        index_char <- index_char + 1
      }
      pos_parenthese_fermee <- index_char-1
      diff_order <- substring(f_x[i_text],pos_newdiff,pos_newdiff)
      f_x[i_text]<-gsub(gsub("(\\[)","\\\\\\1",gsub("(\\])","\\\\\\1",gsub("(\\))","\\\\\\1",gsub("(\\()","\\\\\\1",gsub("(\\+)","\\\\\\1",gsub("(\\*)","\\\\\\1",
                                                                                                                                                substring(f_x[i_text],pos_newdiff-8,pos_parenthese_fermee))))))),
                        paste("(","(",
                              substring(f_x[i_text],pos_newdiff+2,pos_parenthese_fermee),
                              "-(",
                              gsub("t_data\\[t",paste("t_data[t-",diff_order,sep=""),substring(f_x[i_text],pos_newdiff+2,pos_parenthese_fermee-1)),
                              ")",")",
                              sep=""),
                        f_x[i_text])
    }
  }

  #remplacer dataframe[] par submat de RccpArmadillo
  f_x <- as.character(vapply(f_x,function(x)gsub("t_data\\[(t-?\\d?-?\\d?)(,\\d+)\\]","M\\(\\1\\2\\)",x), character(1)))
  f_x <- as.character(vapply(f_x,function(x)gsub("t_data\\[(t-?\\d?-?\\d?-\\(M\\(t-?\\d?,\\d+\\)\\+?\\d?\\))(,\\d+)\\]","M\\(\\1\\2\\)",x), character(1)))

  f_x <- as.character(vapply(f_x,function(x)gsub("f_x_n\\[(\\d+)\\]","f_x_n\\(\\1-1\\)",x), character(1)))


  #remplacer <-  par =
  f_x <- as.character(vapply(f_x,function(x)gsub("<-","=",x), character(1)))

  f_x_f <- c(firstTreated, f_x, lastTreated, "")

  for (i_text in 1:length(f_x_f)){
    while(grepl("\\^\\d",f_x_f[i_text])){ #s'il y a une puissance
      pos_parenthese_ouverte <- integer()
      pos_power <- gregexpr("\\^\\d",f_x_f[i_text])[[1]][1]-2# on se positionne avant de fermer la parenthese de power
      #tant qu'il reste du newdiff a traiter
      #Pour chaque newdiff on veut reperer la parenthese fermante
      nb_parenthese_ouverte <- -1
      index_char <- pos_power
      while (nb_parenthese_ouverte < 0){
        nb_parenthese_ouverte <- nb_parenthese_ouverte +
          (substring(f_x_f[i_text],index_char,index_char)=="(") -
          (substring(f_x_f[i_text],index_char,index_char)==")")
        index_char <- index_char - 1
      }
      pos_parenthese_fermee <- index_char+1
      if(substring(f_x_f[i_text],pos_parenthese_fermee-1,pos_parenthese_fermee-1)=="M"){
        pos_parenthese_fermee=pos_parenthese_fermee-1
      }
      power <- substring(f_x_f[i_text],pos_power+3,pos_power+3)
      f_x_f[i_text]<-gsub(gsub("(\\[)","\\\\\\1",gsub("(\\])","\\\\\\1",gsub("(\\))","\\\\\\1",gsub("(\\()","\\\\\\1",gsub("(\\+)","\\\\\\1",gsub("(\\*)","\\\\\\1",gsub("(\\^\\d+)","\\\\\\1",
                                                                                                                                                                         substring(f_x_f[i_text],pos_parenthese_fermee,pos_power+3)))))))),
                          paste("pow(",
                                substring(f_x_f[i_text],pos_parenthese_fermee,pos_power+1),
                                ",",power,")",
                                sep=""),
                          f_x_f[i_text])
    }
  }
  if(!is.na(cpp_file)){
    if(file.exists(paste(gsub("<.+$","",first[1]),".cpp",sep=""))){file.remove(paste(gsub("<.+$","",first[1]),".cpp",sep=""))}
    cat(f_x_f,file=cpp_file,sep="\n",append = TRUE)
  }else{
    if(file.exists(paste(gsub("<.+$","",first[1]),".cpp",sep=""))){file.remove(paste(gsub("<.+$","",first[1]),".cpp",sep=""))}
    cat(f_x_f,file=paste(gsub("<.+$","",first[1]),".cpp",sep=""),sep="\n",append = FALSE)
  }



}


#' Generate the rcpp file
#'
#' @param model_name name of the model
#' @param p_h_e_bool vector of 3 booleans for prologue, heart, epilogue
#' @param p_h_e_jac list of 3 jacobians for prologue, heart, epilogue
#' @param all_model_vars all variables that are in the model, in alphabetical order
#' @param rcpp_path path where to store the rcpp file
#' @author Charlotte Nudelmann (DG Tresor)
#' @keywords internal
#' @return rcpp source file
create_model_rcpp_source <- function(model_name,p_h_e_bool,p_h_e_jac ,all_model_vars ,rcpp_path){
  rcpp_path<-normalizePath(rcpp_path)
  assertthat::is.dir(rcpp_path)
  ### initialising
  all_model_vars <- all_model_vars[order(all_model_vars)]
  rcpp_file_name <- file.path(rcpp_path,paste0(model_name,"_pomme_newton.cpp"))
  if(file.exists(rcpp_file_name)){file.remove(rcpp_file_name)}

  Rcpp_header <- paste(
'#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends("RcppArmadillo")]]
# include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;
')
  Rcpp_export <- paste('// [[Rcpp::export]]')
  Rcpp_Solver <-paste('// [[Rcpp::export]]
arma::mat Rcpp_solver(arma::mat& M,
                            int& first_date, int& last_date,
                            double& convergenceCriteria,
                            arma::uvec& endo, arma::uvec& prologue_endo, arma::uvec& heart_endo, arma::uvec& epilogue_endo,
                            bool& has_prologue, bool& has_heart, bool& has_epilogue, StringVector& time_line){


  // Delaration des variables
  arma::uvec tt;
  arma::mat Jacobian_n;
  arma::sp_mat sp_jac;
  arma::mat f_x_n;
  arma::vec x_n;
  arma::vec x_n1;
  int t;
  int n;
  double convergence;
  bool s;

  // A chaque periode on fait un newton pour chaque sous partie du modele
  for (t = first_date; t<=last_date ; t++) {
  Rcout << t+1 ;

    // Initialisation a chaque periode
    tt = t;
    M.submat(tt,endo) = M.submat(tt-1, endo);

    //Prologue
    if(has_prologue){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,prologue_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);
      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, prologue_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = prologue_jacobian_f(M, t);
        sp_jac = sp_mat(Jacobian_n);
        f_x_n = prologue_equations_f(M,t);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, prologue_endo) = x_n1.t();
    }

    //Heart
    if(has_heart){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,heart_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);

      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, heart_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = heart_jacobian_f(M,t);
        sp_jac = sp_mat(Jacobian_n);
        f_x_n = heart_equations_f(M,t);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, heart_endo) = x_n1.t();
    }

    //Epilogue
    if(has_epilogue){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,epilogue_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);

      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, epilogue_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = epilogue_jacobian_f(M,t);
        sp_jac = sp_mat(Jacobian_n);
        f_x_n = epilogue_equations_f(M,t);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, epilogue_endo) = x_n1.t();
    }

  }
  return M;
}
')


  solver_head<-'// [[Rcpp::export]]
arma::mat Rcpp_solver(arma::mat& M,
                            int& first_date, int& last_date,
                            double& convergenceCriteria,
                            arma::uvec& endo, arma::uvec& prologue_endo, arma::uvec& heart_endo, arma::uvec& epilogue_endo,
                            bool& has_prologue, bool& has_heart, bool& has_epilogue, StringVector& time_line){


  // Delaration des variables
  arma::uvec tt;
  mat Jacobian_n;
  arma::mat f_x_n;
  arma::vec x_n;
  arma::sp_mat sp_jac;
  arma::vec x_n1;
  int t;
  int n;
  double convergence;
  bool s;

  // A chaque periode on fait un newton pour chaque sous partie du modele
  for (t = first_date; t<=last_date ; t++) {


    // Initialisation a chaque periode
    tt = t;
    M.submat(tt,endo) = M.submat(tt-1, endo);
'
solver_prologue <- '
//Prologue
    if(has_prologue){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,prologue_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);
      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, prologue_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = prologue_jacobian_f(M, t);
        f_x_n = prologue_equations_f(M,t);
        sp_jac = sp_mat(Jacobian_n);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, prologue_endo) = x_n1.t();
    }

'
solver_heart<- '
//Heart
    if(has_heart){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,heart_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);

      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, heart_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = heart_jacobian_f(M,t);
        f_x_n = heart_equations_f(M,t);
        sp_jac = sp_mat(Jacobian_n);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, heart_endo) = x_n1.t();
    }
'

solver_epilogue <- '
//Epilogue
    if(has_epilogue){
      n = 0;
      convergence=convergenceCriteria;
      s = false;
      x_n = M.submat(tt,epilogue_endo).t();
      x_n1 = arma::vec(x_n.size(), arma::fill::ones);

      while (s == false){
        if (n>0) {
          x_n = x_n1;
          M.submat(tt, epilogue_endo) = x_n1.t();
          if(n>10){convergence=10*convergenceCriteria;n=1;}
        }
        Jacobian_n = epilogue_jacobian_f(M,t);
        f_x_n = epilogue_equations_f(M,t);
        sp_jac = sp_mat(Jacobian_n);
        x_n1 = x_n - spsolve(sp_jac, f_x_n);
        s = arma::approx_equal(x_n1,x_n, "absdiff", convergence);
        n++;
      }
      M.submat(tt, epilogue_endo) = x_n1.t();
    }
'
solver_end <- '
 Rcout << "  " ;
 Rcout << time_line(t) ;
 Rcout << " ...  " ;
 }
return M;
}
'
  ##creating source file
  cat(Rcpp_header,file=rcpp_file_name,sep="\n",append = TRUE)

  if((p_h_e_bool)[1] == TRUE){
  cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
  f_x_Rcpp(file.path("temp_paprfn","p_e_f.R"),
           all_model_vars, rcpp_file_name)

  cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
  Jacobian_Rcpp(file.path("temp_paprfn","p_j_f.R"),
           all_model_vars, p_h_e_jac[[1]], rcpp_file_name)
  }

  if((p_h_e_bool)[2] == TRUE){
    cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
    f_x_Rcpp(file.path("temp_paprfn","h_e_f.R"),
             all_model_vars, rcpp_file_name)

    cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
    Jacobian_Rcpp(file.path("temp_paprfn","h_j_f.R"),
                  all_model_vars, p_h_e_jac[[2]], rcpp_file_name)
  }

  if((p_h_e_bool)[3] == TRUE){
    cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
    f_x_Rcpp(file.path("temp_paprfn","e_e_f.R"),
             all_model_vars, rcpp_file_name)

    cat(Rcpp_export,file=rcpp_file_name,sep="\n",append = TRUE)
    Jacobian_Rcpp(file.path("temp_paprfn","e_j_f.R"),
                  all_model_vars, p_h_e_jac[[3]], rcpp_file_name)
  }

 cat(solver_head ,file=rcpp_file_name,sep="\n",append = TRUE)
 if(p_h_e_bool[1]){cat(solver_prologue,file=rcpp_file_name,sep="\n",append = TRUE)}
 if(p_h_e_bool[2]){cat(solver_heart,file=rcpp_file_name,sep="\n",append = TRUE)}
 if(p_h_e_bool[3]){cat(solver_epilogue,file=rcpp_file_name,sep="\n",append = TRUE)}
 cat(solver_end,file=rcpp_file_name,sep="\n",append = TRUE)
 # cat(Rcpp_Solver,file=rcpp_file_name,sep="\n",append = TRUE)

 cat(paste("\n Successfully created the rcpp source files. \n"))
}


