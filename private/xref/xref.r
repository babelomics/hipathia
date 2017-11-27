
### xref functions

process_xref <- function(xfile){
  foreign <- read.table(xfile,header=T,sep="\t",stringsAsFactors=F)
  this_entrezs <- unique(foreign[,1])
  foreign <- foreign[which(foreign[,2]!=""),]
  foreign <- unique(foreign)
  xref <- by(foreign,foreign[,2],function(x) x[,1])
  s <- summary(sapply(xref,length))
  tt <- table(sapply(xref,length))
  return(list(xref=xref,entrezs=this_entrezs,n=nrow(foreign),name=xfile,stats=cbind(do.call("cbind",as.list(s)),n=nrow(foreign)),freqs=tt))
}

process_species <- function(xref_folder){
  xref_files <- list.files(xref_folder,full.names = T)
  xref_files <- xref_files[grep("xref.rdata",xref_files,invert = T)]
  xref <- list()
  all_entrezs <- c()
  stats <- c()
  tts <- list()
  for(xfile in xref_files){
    cat("Processing",xfile,"\n")
    this_xref <- process_xref(xfile)
    print(this_xref$stats)
    xref <- c(xref,this_xref$xref)
    all_entrezs <- c(all_entrezs,this_xref$entrezs)  
    stats <- rbind(stats,this_xref$stats)
    tts[[xfile]] <- this_xref$freqs
  }
  all_entrezs <- unique(all_entrezs)
  entrezs_list <- as.list(all_entrezs)
  names(entrezs_list) <- all_entrezs
  xref <- c(xref,entrezs_list)
  rownames(stats) <- xref_files
  attr(xref,"stats") <- stats
  attr(xref,"freqs") <- tts
  save(xref,file=paste0(xref_folder,"/xref_",xref_folder,".rdata"))
}

### process species

species <- c("hsa","mmu","rno")
for(spe in species){
  cat(">>>Processing",spe,"\n")
  process_species(spe)
  cat("\n\n")
}
