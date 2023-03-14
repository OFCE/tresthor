
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
  Rcpp_Solver <- paste('// [[Rcpp::export]]
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
