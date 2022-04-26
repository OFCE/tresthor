##Help functions


#' Load the Opale network
#'
#' @return opens Opale Network of variables
#' @export
#'
Opale_Net<- function(){
  browseURL(system.file("Opale","opalenet.html",package = "tresthor"))
}

#' Opens the UK example guide
#'
#' @return Opens the UK example guide
#' @export
#'
UK_example<-function(){
  browseURL(system.file("UK_example","guide_exemple_thor.html",package = "tresthor"))
}

#' Opens the Opale in tresthor guide
#'
#' @return Opens the Opale in tresthor guide
#' @export
#'
Opale_guide<-function(){
  browseURL(system.file("Opale","vignetteOpale.html",package = "tresthor"))
}


#' tresthor user manual
#'
#' @return loads the user guide
#' @export
#'
#' @examples
#' tresthor_help()
tresthor_help<-function(){
  browseURL(system.file("vignetteUserManual.html",package = "tresthor"))
}
