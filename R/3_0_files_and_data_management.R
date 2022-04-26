#' Save a thoR model
#' The model is saved as a .rds file.
#' @param model thor.model object to save.
#' @param folder_path directory where to store the .rds. Default : working directory
#'
#' @return rds file stored in location
#' @export
#' @importFrom assertthat is.dir
#' @examples
#' \dontrun{
#' save_model(forecast_model)
#' }
save_model<-function(model,folder_path=getwd()){
  assertthat::is.dir(folder_path)
  folder_path<-normalizePath(folder_path)
  saveRDS(model, file=file.path(folder_path,paste0(model@model_name,".rds",sep=""))  )
}



#' Load a saved model
#'
#' @param model string containing model name (without any extension)
#' @param folder_path directory where to find the model
#' @param file Alternatively provide the complete file path, with the extension. this argument will take over model + folder path argument.
#' @param env environment where to load the model. Default : globalenv()
#' @param source_rcpp Boolean. TRUE if sourcing of the rcpp file is needed
#'
#' @importFrom Rcpp sourceCpp
#' @return loaded model in the environment, ready to be solved.
#' @export
#'
load_model<-function(model=NULL,folder_path=getwd(),file=NULL,env=globalenv(),source_rcpp = TRUE){
  if (!is.null(model)){assertthat::is.string(model)}
  use_file = TRUE
  if(is.null(file) & is.null(model)){stop("Please provide either the full file path, or the name of the model with the folder path.")}
  if(is.null(folder_path)){folder_path <- getwd()}

  if(!is.null(file) & !is.null(model)){cat("Pease provide either file or model and folder_path, not both.\n
                                            The model and folder_path arguments will be ignored.")}
  if(is.null(file) & !is.null(model)){use_file = FALSE}

  if(use_file == FALSE){
    folder_path<-normalizePath(folder_path)
  file <- file.path(folder_path,paste0(model,".rds"))
  }

  assertthat::assert_that(file.exists(file))
  plop<- readRDS(file)
  assertthat::assert_that(class(readRDS(file))=="thoR.model")

  assign(plop@model_name, plop, envir = env )

  if(source_rcpp==TRUE & plop@rcpp ==TRUE){
    cat("\n Loading the rcpp source file... \n")
    Rcpp::sourceCpp(plop@rcpp_source)
  }


}

#' Export the model as a .txt file
#'
#' @param model thoR.model obect
#' @param filename complete path and name of file to save
#'
#' @return Saved file
#' @export
#'
export_model<- function(model, filename ="model.txt"){

  if(class(model)!="thoR.model"){stop("The model must be a thoR.model object.")}
  assertthat::is.dir(dirname(filename))
  equations_list <- model@equation_list
  equations_to_write<-ifelse(equations_list$name==equations_list$id,equations_list$equation,  paste0(equations_list$name,":",equations_list$equation))

  if(grepl(pattern = "\\..*$",filename)==FALSE){filename<-paste0(filename,".txt")}

  if(file.exists(filename)){file.remove(filename)
     cat(paste0("\n Overwriting ",filename))}

  cat("endogenous variables : \n", sep = " ",file = filename,append = TRUE)
  cat(model@endo_list, sep = ",",file = filename,append = TRUE)
  cat('\n ##############', sep = "\n",file = filename,append = TRUE )
  cat('exogenous variables : \n', sep = " ",file = filename,append = TRUE )
  cat(model@exo_list, sep = ",",file = filename,append = TRUE )
  cat('\n ##############', sep = "\n",file = filename,append = TRUE )
  cat('coefficients : \n', sep = " ",file = filename,append = TRUE )
  cat(model@coeff_list, sep = ",",file = filename,append = TRUE )
  cat('\n ##############', sep = "\n",file = filename,append = TRUE )
  cat('equations : ', sep = "\n",file = filename,append = TRUE )
  cat(equations_to_write, sep = "\n",file = filename,append = TRUE )

  cat(paste0("\n The model was saved as ",filename))
}

###Load a sample model file
#' load a sample model file
#'
#' @param path_only TRUE to disable browsing and only show the path instead. Default: FALSE
#'
#' @return Opens a sample model file or the path only if path_only = TRUE.
#' @export
#'
#' @examples
#' \dontrun{model_source_example()}
#' model_source_example(TRUE)
model_source_example<- function(path_only = FALSE){
  if(path_only == FALSE){browseURL(system.file("models", "model_type.txt", package = "tresthor"))
    }else
    {system.file("models", "model_type.txt", package = "tresthor")}
}
