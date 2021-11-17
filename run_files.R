for(ii in 4:length(func_groups)){
  gr <- paste0("G",ii)
  gg <- readLines("G3.rmd")
  
  gg[2] <- gsub("IGHV1-24",names(func_groups[ii]), gg[2])
  
  gg[72] <- gsub("IGHV1-24G3",func_groups[ii], gg[72])
  
  writeLines(gg, con = paste0(gr,".rmd"), sep = "\n")
}

L <- c()
for(ii in 4:length(func_groups)){
  gr <- paste0("G",ii)
  gene <- names(func_groups)[ii]
  
  text <- c('- text: "IGHV1-2"','href: G2.html')
  text[1] <- gsub("IGHV1-2",gene,text[1])
  text[2] <- gsub("G2",gr,text[2])
  
  L <- c(L,text)
}

cat(paste0(L, collapse = "\n"),sep = "\n")

files <- list.files(".","rmd")
for(f in files){
  print(f)
  knitr::knit(f, quiet = T)
}

