## code to prepare `functions_supported` dataset goes here


#
#
# formula_computes <- function(funk){
#   x=sample(c(1:20),size = 1)
#   plus = paste0(funk,"(",x,")")
#   minus = paste0(funk,"(-",x,")")
#   zero = paste0(funk,"(",0,")")
#   raw = paste0(funk,"(x)")
#
#  cat(paste(plus,"=" , eval(parse(text=plus)),"\n",
#            minus,"=" , eval(parse(text=minus)),"\n",
#            zero,"=" , eval(parse(text=zero)),"\n",
#            "derivative \n",raw,":" , Deriv(raw),"\n"
#            ) )
#   }
#
# drule_list <-c("expm1","digamma",  "log2",
#                "log", "log1p", "exp", "logb",
#                 "acosh", "asinh", "sinpi",
#                "cospi",   "asin", "sinh", "cosh",
#                "log10", "atanh", "tanpi",  "acos",
#                "sqrt",  "sin", "atan", "tanh",
#                "abs", "cos",  "tan", "sign")
# drule_list %>% map(~formula_computes(.x))


thor_functions_supported <- c("expm1", "log2",
                              "log", "log1p", "exp", "logb",
                              "acosh", "asinh", "sinpi",
                              "cospi",   "asin", "sinh", "cosh",
                              "log10", "atanh", "tanpi",  "acos",
                              "sqrt",  "sin", "atan", "tanh",
                              "abs", "cos",  "tan", "sign")

usethis::use_data(thor_functions_supported, overwrite = TRUE)



