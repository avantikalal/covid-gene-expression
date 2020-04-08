library(data.table)
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(foreach)
library(doParallel)
library(rtracklayer)

# Function to read GFFs
readGFF=function(gff_file, skipfirst=T){
  gff = import(gff_file)
  if(skipfirst){
    gff = gff[2:length(gff)]
  }
  return(gff)
}


# Function to read fasta > string
readFasta=function(fasta_file, toRNA=F){
  fa = read.fasta(fasta_file, as.string = T, forceDNAtolower = F)
  seq = fa[[1]][[1]]
  if(toRNA){
    seqString = RNAString(gsub("T", "U", seq))
  } else{
    seqString = DNAString(seq)
  }
  seqString = list(seqString)
  names(seqString) = names(fa)
  return(seqString)
}


# Function to read PWMs
readPWMsFromFasta = function (pwm_file) {

  # Read all lines from PWM file
  lines = readLines(pwm_file)

  # Find header lines, start and end of each PWM
  ind = which(substr(lines, 1L, 1L) == ">")
  nseq = length(ind)
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))

  # Get PWM IDs
  ids = lapply(seq_len(nseq), function(i) {
    firstword <- strsplit(lines[ind[i]], "\t")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })

  # Split PWMs
  pwms = lapply(seq_len(nseq), function(i) strsplit(lines[start[i]:end[i]], "\t"))

  # Format as numeric matrix
  pwms = lapply(pwms, function(x) matrix(as.numeric(unlist(x)), ncol=4, byrow=T))

  # Convert to PWMatrix class
  pwms = lapply(seq_len(nseq), function(i){
    PWMatrix(profileMatrix=matrix(c(pwms[[i]]), byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))), ID=ids[[i]])
  })

  # Name with PWM ID
  names(pwms) = ids
  return(pwms)
}

# Multithreaded function to scan sequence(s) with multiple PWMs
ScanSeqWithPWMs = function(seqString, pwmList, seqName, strand="*"){
  sites = foreach(i = 1:length(pwmList), .combine=rbind) %dopar% {
    # Read PWM ID
    id = as.character(names(pwmList)[i])
    # Scan genome
    curr_sites = searchSeq(pwmList[[i]], seqString, min.score="95%", strand=strand, seqname = seqName)
    if(length(curr_sites) > 0){
      # Convert to data table
      curr_sites = as.data.table(writeGFF3(curr_sites))
      if(length(curr_sites) > 0){
        curr_sites[, seq:= curr_sites[, tstrsplit(attributes, split=";|=", perl=T)][, V6]]
        curr_sites[, attributes:=NULL]
        curr_sites[, Matrix_id:= id]
      }
    }
  }
  return(sites)
}

# Simulation function

simulateSeq = function(freqs, N){
  bases=c()
  for(base_type in names(freqs)){
    bases = c(bases, rep(base_type, freqs[base_type]))
  }
  sim_genomes = list()
  for(i in 1:N){
    sim_genomes[[i]] = DNAString(paste(sample(bases), collapse=""))
  }
  sim_genomes = DNAStringSet(sim_genomes)
  names(sim_genomes) = paste0("sim", 1:N)
  return(sim_genomes)
}