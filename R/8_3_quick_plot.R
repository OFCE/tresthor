##quick_plot
#' Draw a line plot quickly
#' @param variables variables to plot
#' @param database data.frame
#' @param start first date to plot Default: NULL (first observation)
#' @param end last date to plot. Default: NULL (last observation)
#' @param index_time Name of the column which contains the time info. Default: "date"
#' @param growth_rate TRUE if growth rates should be plotted instead of levels. Default : FALSE
#' @param title Title of the graph
#' @param colours vector of colours for the plotted line, must be at least as many colours as variables. extra colours will be ignored. If not enough colours are provided, the colours will revert to the defaul ones.
#' @import scales
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @return ggplot graph
#' @export
#'
quick_plot<-function(variables,database=data,start=NULL,end=NULL,index_time="date",growth_rate=FALSE,title=NULL,colours=NULL){
  options(scipen = 99)
  assertthat::assert_that(inherits(database,"data.frame"))
  if(!is.null(title)){assertthat::is.string(title)}else{title <- ""}
  if(!index_time %in% colnames(database)){stop('index_time must be a variable of the database.')}

  if(!is.null(start)){if(!start %in% database[,index_time]){stop("start must be an element of index_time.")}}
  if(!is.null(end)){if(!end %in% database[,index_time]){stop("end must be an element of index_time.")}}

  if(is.null(start)){start <- head(database[,index_time],1)}
  if(is.null(end)){end <- tail(database[,index_time],1)}
  if(start > end){stop("start must be before end.")}

  if(prod(variables %in% colnames(database) )==0){stop("Some variables specified could not be found.")}
  database$plop <- c(1:nrow(database))
  graph_data <- database %>% select(all_of(index_time),all_of(variables),plop)

  istart <- which(database[,index_time]==start)
  iend <- which(database[,index_time]==end)
  if(growth_rate){
    database <- database %>% mutate_at(variables,.funs = ~(.x/lag(.x) -1))
  }
  graph_data <- database %>% filter(plop >= istart & plop <= iend)
  if(inherits(graph_data[,index_time],"Date")){
    graph_data[,index_time]<-reformat_date(graph_data[,index_time])
  }

  graph_data$plopsy<-graph_data[,index_time]

  data_courbes <- graph_data %>% tidyr::pivot_longer(all_of(variables),"variable" )
  n_obs<-nrow(graph_data)

  if(n_obs>=14){
    angle_date <- 90
  }else{
    angle_date <- 0}

  n_vars<- length(variables)

  if(!is.null(colours)){
    if(length(colours) >= n_vars){
      colour_vec <- head(colours, n_vars)
    }else{
      thor_palette<- c("blue3","darkorange","turquoise3",
                       "palevioletred3","gold2","firebrick3",
                       "springgreen4","deeppink4","sienna4","yellowgreen",
                       "peachpuff3","indianred3")

      if (n_vars <= length(thor_palette)) {
        colour_vec <- head(thor_palette, n_vars)
      } else{
        colour_vec <- palette(rainbow(n_vars))
      }

    }
  }else{
    thor_palette<- c("blue3","darkorange","turquoise3",
                     "palevioletred3","gold2","firebrick3",
                     "springgreen4","deeppink4","sienna4","yellowgreen",
                     "peachpuff3","indianred3")

    if (n_vars <= length(thor_palette)) {
      colour_vec <- head(thor_palette, n_vars)
    } else{
      colour_vec <- palette(rainbow(n_vars))
    }

  }


  if (growth_rate) {
    if (max(na.omit(data_courbes$value)) > 0.1) {
      break_seq <- seq(-100, 100, 0.02)
    }else{
      if (max(na.omit(data_courbes$value)) > 0.015) {
        break_seq <- seq(-100, 100, 0.01)
      }else{
        break_seq <- seq(-100, 100, 0.0025)
      }
    }
  }

  if(growth_rate){
    plot<-ggplot(data=data_courbes,aes(group= variable,y = value, x = plopsy,color= variable)) +
      geom_point(size=0)+
      geom_line()+

      scale_color_manual(values=colour_vec,name = "")+

      scale_y_continuous(name = "",labels = scales::percent_format(accuracy = 0.25),breaks=break_seq, minor_breaks = break_seq )+
      geom_hline(yintercept = 0,color = "grey16")+
      ggtitle(paste(title))+

      labs( x = " ")+
      theme(legend.position="bottom",axis.text.x = element_text(angle = angle_date) ,panel.border = element_rect(colour = "black", fill=NA, size=0.7),panel.background = element_rect(fill = 'white', colour = 'grey78'),panel.grid = element_line(color="grey89")) +
      guides(color = guide_legend(ncol = 2, order = 1))

  }else{
    plot<-ggplot(data=data_courbes,aes(group= variable,y = value, x = plopsy,color= variable)) +
      geom_point(size=0)+
      geom_line()+
      scale_color_manual(values=colour_vec,name = "")+

      geom_hline(yintercept = 0,color = "grey16")+
      ggtitle(paste(title))+

      labs(y= " ", x = " ")+
      theme(legend.position="bottom",axis.text.x = element_text(angle = angle_date) ,panel.border = element_rect(colour = "black", fill=NA, size=0.7),panel.background = element_rect(fill = 'white', colour = 'grey78'),panel.grid = element_line(color="grey89")) +
      guides(color = guide_legend(ncol = 2, order = 1))
  }

  plot
}
