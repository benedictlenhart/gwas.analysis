.libPaths(c("/project/berglandlab/Rlibs_4.3.1/")); .libPaths()
#.libPaths(c("/home/bal7cg/R/goolf/4.3")); .libPaths()
library(RhpcBLASctl)
blas_set_num_threads(1)
#***************************
### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(5)

#library(ggplot2)

print("libraries loaded")

args = commandArgs(trailingOnly=TRUE)
#create job id from array
id = as.numeric(args[1])
#id = 1
#id = 1#use this when testing to see if script functions correctly
numCores <- as.numeric(args[2])
#numCores = 4
print(id)
print(numCores)

# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/standard/vol186/bergland-lab/Adam/gwas/"#path of working directory

print(INPUT)
#load in reftable
reftable = readRDS(paste0(INPUT,"dgrpreftable"))

#load in array table
arrayref = readRDS(paste0(INPUT, "ref.table"))
#filter to this array's section
arrayref = arrayref[array == id]
print(arrayref)
out = foreach(f = as.vector(arrayref[[1]]), .errorhandling = "remove") %do% {
  # f = 1
  print(f)
  #load in our permuted phenotype file based on array
  phenotable = readRDS(paste0(INPUT,"/phenoperms/phenopermutation",f))
  #melt phenotable
  melt = melt(phenotable, id.vars = ("DGRP"), variable.name = "fullpheno", value.name = "mean")
  #phenotype determined by array table
  pheno.id = as.character(arrayref[ID_var_1 == f]$ID_var_2)


  dt = melt[fullpheno == pheno.id]
  #merge in the values for reftable
  dt = merge(dt, reftable, by = "DGRP" )
  dt = na.omit(dt)#remove na's
  #fix ral ids
  dt = dt %>%
    dplyr::rename("ral_id" = DGRP) %>%
    mutate(ral_id = gsub("DGRP", "line", ral_id))
  print(dim(dt))





  #re organize
  colnames(dt)[3] <- "avg"




  #try making to numberic
  dt$ral_id = gsub( "line_","", dt$ral_id)
  dt$ral_id2 = as.numeric(dt$ral_id)#numeric removes excess zeroes
  dt$ral_id2 = paste0("line_", dt$ral_id2)

  #create a second loop that goes through chromosomes
  chroms= c("2L", "2R", "3L", "3R", "X")
  chrom.out = foreach( c = chroms) %do% {
    #c = "2L"
    #load in the grm for that chromosome
    grm = readRDS(paste0(INPUT, "No", c, "GRM"))
    #assign gds filename based on chromosome
    ingds = paste0(INPUT,c,".gds")
    modelqtl <- glmmkin(fixed = avg ~ WolbachiaStatus_NA, data = dt, kins = grm, id = "ral_id2",family = gaussian(link = "identity"))

    if(!file.exists(paste(INPUT,"permlocoGWAS/perm", f, sep=""))) {
      message("making directory")
      dir.create(paste(INPUT,"permlocoGWAS/perm", f, sep=""))
    } else {
      message("directory exists")
    }

    outputs = paste(INPUT ,
                    "permlocoGWAS/perm", f, "/",
                    f,"-",pheno.id,"_",c, ".txt",
                    sep = "")
    glmm.score(modelqtl, infile = ingds, outfile = outputs, MAF.range = c(0.05,0.95), miss.cutoff = 0.15, ncores = numCores ,
               nperbatch = 100, verbose = T)

    print ("done")
    # #replace file with rds file to save space
    x = fread(outputs)
    saveRDS(x, paste(outputs,".RDS", sep=""))
    file.remove(outputs)#remove txt file
    # print("finished")

  }


}
