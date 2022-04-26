mainf<-function(ar,ma){
  n_ar<-length(ar)
  n_ma<-length(ma)

  sor<-vector()
  sor[1]<-ma[1]  ## le premier terme est simplement egal au premier terme de la partie ma */
  i<-2  ## on traite donc maintenant a partir du second terme */
  while( i<=n_ma ) ## les termes ma n'interviennent directement dans le calcul du coefficient que tant qu'ils existent */
  {
    nbmin<-min(i-1,n_ar)
    somme<-0
    for (j in 1:nbmin)
    {
      somme<-somme+sor[i-j]*ar[j]
    }
    sor<-c(sor,ma[i]+somme)
    i <- i +1
  }
  for (i in n_ma+1:100)
  {
    nbmin<-min(i-1,n_ar)
    somme<-0
    for (j in 1:nbmin)
    {
      somme<-somme+sor[i-j]*ar[j]
    }
    sor<-c(sor,somme)
  }


  return(sor)

}

contrib<-function(serie,serie_mainf){
  longueur_originale <- length(serie)
  ###verifie que la serie ne finisse pas en NA
  if (is.na(serie[length(serie)])==TRUE){stop("series ends with NAs.")}

  ##prend la serie sans NA
  serie2<-na.omit(serie)
  debut <- 1
  fin<-length(serie2)  ##date de fin de la serie */
  n_mainf<-length(serie_mainf)
  sortie<-vector()  ##initialisation a  blanc de sortie */
  nb_NA <- longueur_originale-fin
  for (i in debut:fin)
    ##bouclage jusqu'a la date de fin de la serie en entree */
  {
    nbret0<-min(i,n_mainf)
    deb<-max(1,(i-n_mainf+1))
    ## nombre de retards sur lequel peut être calculee la contrib, compte tenu de la date de depart de la serie et du nombre
    ## de la date i sur laquelle est calculee cette contrib

    serie_mainf_souspart<-rev(serie_mainf[1:nbret0])
    ## extraction de la partie du mainf dont on a besoin et inversion

    serie_souspart<-serie2[deb:i]

    ## extraction de la partie de la serie dont on a besoin

    serie_calcul<-(serie_souspart %*% serie_mainf_souspart)[1,1]
    ## multiplication de la structure de mainf par la structure de retard de la serie

    sortie<-c(sortie,serie_calcul) ##ajout aux calculs precedents */
  }
  contrib_serie<- c(rep(NA,nb_NA),sortie)##refomatage de la serie pour la ramener a la date voulue */
  return(contrib_serie)
}



derivative_time_table<-function(formula_num_coeff,maxlag){
options(scipen = 999)

f<-gsub("\\\n","",formula_num_coeff)
f<-gsub("\\s+","",f)

f<-gsub("(mylg|lag)\\((\\w+),-?([0-9]+)\\)","lagxxplopxx\\2xxplopxx\\3",f)
f<-gsub("log\\((\\w+(xxplopxx\\w+)*)\\)","logxxplopxx\\1",f)
f<-gsub("delta\\(([0-9]+),(\\w+(xxplopxx\\w+)*)\\)","(\\2-lagxxplopxx\\2xxplopxx\\1)",f)
f<-gsub("((\\=|\\+|\\-|\\*|/|\\()([0-9]+(\\.[0-9]+)?)(\\=|\\+|\\-|\\*|/|\\)))","MOOGLES\\1MOOGLES",f)
f<-gsub("((MOOGLES)([0-9]+(\\.[0-9]+)?)(\\=|\\+|\\-|\\*|/|\\)))","MOOGLES\\1MOOGLES",f)
f<-gsub("lagxxplopxx(logxxplopxx)?lagxxplopxx(\\w+)xxplopxx([0-9]+)xxplopxx([0-9]+)","lagxxplopxx\\1\\2xxplopxx\\3+\\4",f)

f<-gsub("logxxplopxxlagxxplopxx","lagxxplopxxlogxxplopxx",f)
f<-gsub("lagxxplopxxlagxxplopxx","lagxxplopxx",f)
f<-gsubfn::gsubfn("(xxplopxx\\w+)xxplopxx([0-9]+\\+[0-9]+)", function(n,m) paste0(n, "xxplopxx", eval(parse(text = m))), f)
f<-gsub("MOOGLES","",f)
f<-gsub("=","-(",f)
f<-gsub("$",")",f)
f<-gsub("--","+",f)
f<-gsub("-\\+","-",f)
f<-gsub("\\+-","-",f)
f<-gsub("\\+\\+","+",f)

list_var_deriv<-get_variables_from_string(f)
list_var_deriv<-gsub('xxplopxx',".",list_var_deriv)
deriv_f<-gsub('xxplopxx',".",f)
# if((stringr::str_count(deriv_f,"\\(") < stringr::str_count(deriv_f,"\\)")) & grepl("\\)$",deriv_f)){
#   deriv_f<-gsub("\\)$","",deriv_f)
# }
#### Liste les variables et les versions laggees
T_0<-list_var_deriv
T_0<-gsub("(\\.[0-9]+)$","",T_0)
T_0<-gsub("^lag\\.","",T_0)
T_0<-unique(T_0)

table_deriv<-data.frame(t_0=T_0,row.names = gsub("log\\.","",T_0) )
for ( i in 1:maxlag){
  table_deriv[,paste0("t_",i)]<- paste("lag",table_deriv[,1],i,sep=".")
}


#### Calcul des derivees partielles
derivees_partielles<-table_deriv
for(i in colnames(derivees_partielles)){
  for(j in rownames(derivees_partielles)){
    derivees_partielles[j,i]<-Deriv(deriv_f,table_deriv[j,i])
  }
  derivees_partielles[,i]<-as.numeric(derivees_partielles[,i])
}
res <- cbind(name_log=table_deriv[,1],derivees_partielles)
if(sum(is.na(res))>0){stop("Try to simplify the formula by using one variable instead of expressions of multiple variables.")}else{return<-res}
}

calcul_contrib<-function(var, data,deriv_table,mainf_endo,index_time){
  deriv_table_num <- select(deriv_table,-c("name_log"))

  mainf_exo<- -1 * (deriv_table_num[var,])
  while (mainf_exo[length(mainf_exo)]==0) {mainf_exo[length(mainf_exo)]<-NULL}
  mainf_vec <- mainf(as.matrix(mainf_endo),as.matrix(mainf_exo))

  log.true <- grepl("^log\\.",deriv_table[var,1])
  if (log.true==TRUE){
    data[,paste0(var,".contrib")]  <- contrib(delta(1, log(data[,var])),mainf_vec)
  }else{
    data[,paste0(var,".contrib")]  <- contrib(delta(1, data[,var]),mainf_vec)
  }
  return <- data[,c(index_time,paste0(var,".contrib"))]
}

quarterly_to_annual<-function(series,database,index_year,quarter_start=1,multiply=100){
  if(class(database)!="data.frame") {
    stop('Database must be a data.frame.')
  }
  if (class(series) != "character") {
    stop("serie must be character object of length 1.")
  }
  if (length(series) != 1) {
    stop("Only 1 series at time.")
  }
  if (!series %in% names(database)) {
    stop("Series not found in database.")
  }
  name_series <- paste0(series, ".y")


  if (length(database[, series]) != length(index_year)) {
    stop("The year vector needs to as long as the series vector")
  }
  if (length(database[, series]) < 8) {
    stop("Il faut des series plus longues (minimum 8 observations)")
  }

  if (quarter_start == 1) {
    a <- c(1, 2, 3, 4)
  } else
  {
    if (quarter_start == 2) {
      a <- c(2, 3, 4, 1)
    } else
    {
      if (quarter_start == 3) {
        a <- c(3, 4, 1, 2)
      } else
      {
        if (quarter_start == 4) {
          a <- c(4, 1, 2, 3)
        } else
        {
          stop("quarter start should be 1. 2. 3. or 4.")
        }
      }
    }
  }
  b<-length(database[, series])
  reps<-ceiling(b / 4)
  trim_vect<-c(rep(a, reps))
  trim_vector<-trim_vect[1:b]
  x <-database[, series]
  x0=1+x
  y0=1/x0
  z0=x0*(1 + lead(x0, 1) * (1 + lead(x0, 2) * (1 + lead(x0, 3)))) / (1 +
                                                                       lag(y0, 1) * (1 + lag(y0, 2) * (1 + lag(y0, 3))))

  tann_data<-cbind((z0 - 1) * multiply, index_year, as.numeric(trim_vector))
  rownames(tann_data)<-tann_data[, 2]
  colnames(tann_data)<-c(name_series, "year", "trim_vector")
  tann_data_fin<-as.data.frame(tann_data) %>% dplyr::filter(trim_vector ==
                                                              1) %>% dplyr::select(-trim_vector)

}

