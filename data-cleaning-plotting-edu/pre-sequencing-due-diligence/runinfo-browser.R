library(tidyverse)
library(rio)
setwd("~/projects/bioinfo/echinomics/strongylocentrotus_purpuratus/sequencing/spur5.0")

runinfo = import("runinfo.csv")
names.runinfo = runinfo %>% names()

browse.runinfo = function(var.runinfo){
  
  list.var = runinfo %>% select(var.runinfo) %>% unique()  
  return(list.var)
  
}