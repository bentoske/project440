# this script is written to run a multiple sequence alignment on orthologous DNA data between species and
# that alignment is used to create a phylogenetic tree. Currently, this is a work in progress. 

# this script may also generate a few 502 errors, but will run if you wait them out. 
# particularly during rMSA package install 

# sample dataset is provided in current github repository as BRCA1.fasta

# PROGRESS 
# [x] load multifasta file
# [x] align using MAFFT 
# [x] generate a tree as an output 
# [ ] tree shows species instead of accession numbers 
# [ ] convert script into function 
# [ ] apply script to multiple orthologous species to obtain multiple trees 

# OPTIONAL UPGRADES 
# [ ] run multiple versions of MAFFT for more or less than 200 sequences OR 
# [ ] have multiple options for MSA 
# [ ] run from the shell 
# [ ] generate trees using multiple methods with the ability for user to choose 

# install packages + load libraries - these are sequential. If all installed, just load libraries. 
install.packages("remotes") # allows for install from github
install.packages("BiocManager") # bioconductor
library(remotes) 
library(BiocManager)
remotes::install_github("mhahsler/rMSA") # contains MAFFT
# update all packages when prompted 
# may need to force if changing R versions 
library(rMSA)
install.packages("seqinr")
BiocManager::install("msa") # not using clustal but need msa convert function
install.packages("taxize")
library(ape) # needed for tree
library(seqinr) # needed for tree
library(msa) # needed for tree
library(taxize) # used to convert scientific names 

# install.packages("ggmsa") # unsure if i'm using this 
# library(ggmsa)

# read multi-FASTA
dna <- readDNAStringSet("BRCA1.fasta", format="fasta",
                 nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

# this series of loops will separate the names from accession numbers and other header information. 

names(dna) -> names
count <- length(names) 
for(i in 1:count){
  newname.split <- c()
  oldname.split <- as.vector(t(as.data.frame(strsplit(names[[i]], " "))))
  length <- length(oldname.split)
  w <- 2 # this will account for skipping accession number 
  while(w < 6){
    if(oldname.split[w] == "PREDICTED:"){
      w <- w+1
    }else if(oldname.split[w] == "BRCA1" || oldname.split[w] == "breast" || oldname.split[w] == "BRCA1,"){
      w <- 6
    }else{
      newname.split <- c(newname.split, oldname.split[w])
      w <- w+1
    }
  }
  newname <- paste(newname.split, collapse = '_')
  
  names[i] <- newname
}
# will convert scientific names into common names - this needs revision since it doesn't find all of them
names <- sci2comm(names)
names <- gsub(" ", "_", names)

names(dna) <- names
names
# run alignment 
al <- mafft(dna) # this may take a minute depending on number of sequences

# run to see alignment - optional
detail(al)

# convert multisequence alignment
al <- msaConvert(al, type="seqinr::alignment")

# make tree using neighbor joining method
d <- dist.alignment(al,"identity")
tree <- ladderize(nj(d))

# plots tree 
plot(tree, main = "BRCA1 in 33 Species") # change title later 

# ladderize
tr33 <- ladderize(tree)

plot(tr33, main = "BRCA1 in 33 Species")

