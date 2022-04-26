dictionnaire_variables_opale<-read.csv2("essais/dico_opale_brut") %>%  rename(nom_long = Nom) %>%
  rename_all(~tolower(.x)) %>%
  mutate(categorie = ifelse(categorie=="","Variable",categorie),
         nom_long = ifelse(grepl("coûts salariaux unitaires|prix à l'export|taux de change",nom_long),paste0(nom_long," - ",precision),nom_long ),
         nom_long = ifelse(categorie=="Coefficient économétrique",paste0(nom_long," - ",precision),nom_long),
         nom_long = ifelse(code=="cice_d39","Crédit d'impôt pour la compétitivité et l'emploi (série jusqu'à 2018)",nom_long),
         nom_long = gsub("(^[a-z])","\\U\\1",nom_long,perl = TRUE),
         nom_court = ifelse(categorie=="Variable" & grepl("^td_",code),gsub("\\(prix chaînés\\)|en volume|en volume chaîné| en valeur|\\(prix année N-1\\)" ,"", nom_long) ,nom_long),
         nature = ifelse(categorie=="Variable" & grepl("^td_|^tc_",code)& grepl("(0|1|7|7_ch|3|5_ch)$",code),code,""),
         nature = gsub("^.*0$","prix année N-1",nature ),
         nature = gsub("^td_.*7$","volume",nature ),
         nature = gsub("^tc_.*7$","Equivalent temps plein",nature ),
         nature = gsub("^tc_.*1$","Personnes physiques",nature ),
         nature = gsub("^.*7_ch$","volume chaîné",nature ),
         nature = gsub("^.*5_ch$","prix chaînés",nature ),
         nature = gsub("^.*3$","valeur",nature ),

         nature = ifelse(categorie=="Variable" & grepl("^contpib",code),code,nature),
         categorie = ifelse(categorie=="Variable" & grepl("^contpib",code),"Contribution",categorie),
         nature = gsub("^contpib3.*$","Contribution au PIB valeur",nature ),
         nature = gsub("^contpib7.*$","Contribution au PIB volume",nature ),
         nature = gsub("^contpib5.*$","Contribution au déflateur du PIB",nature )
  )

rownames(dictionnaire_variables_opale)<-dictionnaire_variables_opale$code
usethis::use_data(dictionnaire_variables_opale, overwrite = TRUE)
