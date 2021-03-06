library(plyr)
library(dplyr)   #filter, select
library(magrittr)  #%<>%
library(stringr)
library(purrr)     # %>%
library(functional)

#setwd("C:/Users/Youngjun Park/Desktop/leb/bipolar_STR/")
old.data <- read.table("old_data.txt")
alleles <- read.table("all_chr.converted.length.txt")
alleles.1 <- select(alleles, 1, 2, 3)
alleles.1$V1 <- paste(alleles.1$V1, "_first_chr", sep="")
chr1.row <- alleles.1$V1 %>% unique
chr1.col <- select(alleles, 2) %>% unique %>% reduce(c) %>% sort
#chr1.names <- rep(chr1.names, 2)

alleles.2 <- select(alleles, 1, 2, 4)
alleles.2$V1 <- paste(alleles.2$V1, "_second_chr", sep="")
chr2.row <- alleles.2$V1 %>% unique
chr2.col <- select(alleles, 2) %>% unique %>% reduce(c) %>% sort

add.rownames <- function(missing.array)
{
  rownames(missing.array) <- missing.array[,1]
  return(missing.array)
}

add.colnames <- function(missing.names, missing.array)
{
  missing.array[[missing.names]][1] <- map(missing.array[[missing.names]][1], as.integer) %>% reduce(c)
  colnames(missing.array[[missing.names]]) <- missing.names
  return(missing.array[[missing.names]])
}

alleles.1.split <- split(alleles.1, alleles$V2) 
alleles.1.split <- alleles.1.split[names(alleles.1.split) != "D4S395" & names(alleles.1.split) != "D17S802"] %>% map(add.rownames) %>% map(select, 3) 
alleles.1.split <- map(names(alleles.1.split), add.colnames, alleles.1.split) %>% data.frame

alleles.2.split <- split(alleles.2, alleles$V2)
alleles.2.split <- alleles.2.split[names(alleles.2.split) != "D4S395" & names(alleles.2.split) != "D17S802"] %>% map(add.rownames) %>% map(select, 3) 
alleles.2.split <- map(names(alleles.2.split), add.colnames, alleles.2.split) %>% data.frame

alleles.combined <- rbind(alleles.1.split, alleles.2.split) %>% t
new.data <- alleles.combined[order(rownames(alleles.combined)), order(colnames(alleles.combined))] %>% data.frame
#write.table(new.data, "new_chr1.input")

#new.compare <- t(new.data) %>% apply(2, is.na) %>% apply(2, which) %>% map(length) %>% data.frame %>% t

#filter by comparing with old data
#assign NA to new and old data that have alleles present in one and but not in the other
old.data.filtered <- filter(old.data, is.element(rownames(old.data), rownames(new.data)))
rownames(old.data.filtered) <- rownames(new.data)

new.data[is.na(old.data.filtered)] <- NA

old.data.filtered[is.na(new.data)] <- NA

not.na <- function(array)
{
  return(!is.na(array))
}

new.uniq <- apply(new.data, 1, unique) %>% map(Compose(not.na, length))
old.uniq <- apply(old.data.filtered, 1, unique) %>% map(Compose(not.na, length))

new.test <- new.data[(map2(new.uniq, old.uniq, `==`) %>% reduce(c)),]
old.test <- old.data.filtered[(map2(new.uniq, old.uniq, `==`) %>% reduce(c)),]

#return true if its odd else false
check.even <- function(value)
{
  if(value == -1)
    return(FALSE)
  else if(as.integer(value) %% 2 == 0)
    return(FALSE)
  else
    return(TRUE)
}
#return true if its even else false
check.odd <- function(value)
{
  if(value == -1)
    return(FALSE)
  else if(as.integer(value) %% 2 == 0)
    return(TRUE)
  else
    return(FALSE)
}

filter.by.allele.distance <- function(arr)
{
  num <- arr[1]
  for (x in arr)
  {
    if (x != -1)
      num <- x
      break
  }
  if(num %% 2 == 0)
    val <- map(arr, check.even) %>% unlist %>% which %>% length
  else
    val <- map(arr, check.odd) %>% unlist %>% which %>% length
  if(val == 0)
    return(TRUE)
  else
    return(FALSE)
}

new.test[is.na(new.test)] <- -1
old.test[is.na(old.test)] <- -1

new.df <- new.test[apply(new.test, 1, filter.by.allele.distance),]
old.df <- old.test[apply(old.test, 1, filter.by.allele.distance),]
new.df.filtered <- new.df[is.element(rownames(new.df), rownames(old.df)),]
old.df.filtered <- old.df[is.element(rownames(old.df), rownames(new.df)),]

new.df.filtered[new.df.filtered == -1] <- NA
old.df.filtered[old.df.filtered == -1] <- NA

write.table(new.df.filtered, "new_str.txt")
write.table(old.df.filtered, "old_str.txt")

#is.element(rownames(new.df.filtered), rownames(old.df.filtered))
#is.element(rownames(old.df.filtered), rownames(new.df.filtered))

