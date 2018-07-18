library(plyr)
library(dplyr)   #filter, select
library(magrittr)  #%<>%
library(stringr)
library(stringi)
library(purrr)     # %>%
library(functional)


template <- read.table("new_str.txt")
sample.names <- colnames(template) %>% str_replace("_first_chr", "") %>% str_replace("_second_chr", "") %>% unique

old.str <- read.table("old_str.data")
new.str <- read.table("new_str.data")
rownames(old.str) <- old.str[,1]
old.str <- old.str[,-1]
rownames(new.str) <- new.str[,1]
new.str <- new.str[,-1]

new.str.compare <- map2(new.str[c(T,F)], new.str[c(F,T)], paste, sep = " ") %>% data.frame
#new.str.compare <- map2(new.str[c(T,F)], new.str[c(F,T)], paste, sep = " ") %>% map(map, stri_sort) %>% map(unlist) %>% data.frame
rownames(new.str.compare) <- rownames(new.str)
colnames(new.str.compare) <- sample.names

old.str.compare <- map2(old.str[c(F,T)], old.str[c(T,F)], paste, sep = " ") %>% data.frame
#temp <- old.str.compare %>% map(map, mixedsort) %>% map(unlist) %>% data.frame
rownames(old.str.compare) <- rownames(old.str)
colnames(old.str.compare) <- sample.names

compare.str <- function(arr1, arr2)
{
 return(map2(arr1, arr2, is.element) %>% unlist)
}

find.concordance <- function(arr1, arr2)
{
  match <- map2(arr1, arr2, compare.str) %>% unlist %>% which %>% length
  total <- dim(arr1)[1] * dim(arr1)[2]
  return(match/total)
}

#####################################################
str.concordance <- find.concordance(old.str.compare, new.str.compare)###
str.concordance[2] <- nrow(old.str.compare)
names(str.concordance) <- c("Concordance", "STR")
#####################################################
write(str.concordance, "concordance.out")


