#' Detect frequency
#'
#' @param date_vector A vector of dates (class "Date").
#' @keywords internal
#' @return one of (undetermined,daily, weekly, monthly, quarterly, biannual or annual)
detect_frequency<-function(date_vector){

  if(inherits(date_vector,"Date")==FALSE){stop("The date_vector must be of class 'Date'.")}
  if(length(date_vector)<2){stop("At least 2 dates must be in the vector")}
  if(sum(is.na(date_vector))>0){stop("The date_vector contains some NAs.")}

  diff<-(date_vector - dplyr::lag(date_vector))

  mean_diff<-mean(diff[!is.na(diff)])

  if(mean_diff == 7 ){freq <- "weekly"}else{
    if(mean_diff > 27 & mean_diff < 32){freq <- "monthly"}else{
      if(mean_diff > 88 & mean_diff < 93){freq <- "quarterly"}else{
        if(mean_diff > 179 & mean_diff < 183){freq <- "biannual"}else{
          if(mean_diff > 363 & mean_diff < 368){freq <- "annual"}else{
            if(mean_diff < 5 ){freq <- "daily"}else{
              freq<-"undetermined"
            }
          }
        }
      }
    }
  }

  return<-freq
}


#' Reformat a date vector based on its frequency
#'
#' @param date_vector A vector of dates (class "Date").
#'
#' @return a character vector made to look nicer for graphs
#' @export
#'
#' @examples
#' date<-seq.Date(from=as.Date("2019-02-04"),by = "quarter",length.out = 10)
#' print(reformat_date(date))
reformat_date<-function(date_vector){
  guess_frequency <- detect_frequency(date_vector)

  date_vector_bis<-date_vector

  if(guess_frequency == 'monthly'){
    date_vector_bis<- as.character(format(date_vector, "%b %Y"))}

  if(guess_frequency == 'quarterly'){
    quarters <- quarters(date_vector)

    date_vector_bis<- paste0(format(date_vector, "%Y"),quarters)  }

  if(guess_frequency == 'biannual'){
      month1<- format(date_vector[1], "%m")
      month2<- format(date_vector[2], "%m")
      if (month1<month2){seq_S<- c("S1" , "S2")}else{seq_S<- c("S2" , "S1")}
      date_vector_year<-format(date_vector, "%Y")

      date_vector_bis<-paste0(date_vector_year,seq_S)
  }

  if(guess_frequency == "annual"){
    date_vector_bis<- format(date_vector, "%Y")
  }

  return<-as.character(date_vector_bis)

}
