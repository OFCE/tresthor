#' @title compute n order differences while keeping the vector same length as original one
#' @param n number of periods to diff
#' @param x numeric vector of variable to diff
#' @return numeric vector of diff
#' @rdname newdiff_delta
#' @export
newdiff<-function(n=1,x){
  diff.x<-diff(x,lag=n)
  return<-c(rep(NA,n),diff.x)
}

#' @rdname newdiff_delta
#' @export
delta<-function(n=1,x){
  diff.x<-diff(x,lag=n)
  return<-c(rep(NA,n),diff.x)
}



#' @title Add coefficients to the dataset
#' This function adds coefficients from a table to the dataframe that will be used by the solver
#'
#' @param listcoeff data frame with at least two columns : 1 for name of coeff , another one for value of coeff
#' @param database target data base, Default: t_data
#' @param pos.coeff.name which column in list coeff has the name, Default: 1
#' @param pos.coeff.value which column in list coeff has the value, Default: 2
#' @param overwrite boolean. Whether to overwrite existing values of coefficients already in the database. Default : TRUE.
#' @return returns the original dataset + the coeffs
#' @rdname add_coeffs
#' @export
add_coeffs<-function(listcoeff,database=t_data,pos.coeff.name=1,pos.coeff.value=2,overwrite=TRUE){

  data<-database

  ncoeffs=nrow(listcoeff)
  size_data<-nrow(database)

  listcoeff[,pos.coeff.name]<-tolower(listcoeff[,pos.coeff.name])

  listest<-listcoeff[,pos.coeff.name]
  tdata_list<-colnames(database)

  newcoeffs<-setdiff(listest,tdata_list)
  existingcoeffs<-intersect(listest,tdata_list)

  lco<-as.data.frame(t(matrix(listcoeff[,pos.coeff.value] , ncoeffs , size_data )))
  colnames(lco)<-listcoeff[,pos.coeff.name]

  if (overwrite==TRUE){
    data1<-data[,!(names(data) %in% existingcoeffs)]
    fulldata<-cbind(data1,lco)
  }else{
    if(length(newcoeffs)==0){stop(print="No new coefficients to add.")}
    fulldata<-cbind(data,lco[,newcoeffs])
    colnames(fulldata)<-c(colnames(data),newcoeffs)}

  return <-fulldata
}
