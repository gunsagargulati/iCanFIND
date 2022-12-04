
#Functions
#Function
calculateMutType <- function(sdata, geneName){
  polyphen <- unlist(lapply(sdata, function(x) ifelse(length(which(gsub(".*_", "", x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(geneName)),]$Polyphen_Prediction) %in% c("damaging")))>0, 1, 0)))
  sift <- unlist(lapply(sdata, function(x) ifelse(length(which(x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(geneName)),]$SIFT_Prediction %in% c("deleterious")))>0,1, 0)))
  definite <- unlist(lapply(sdata, function(x) ifelse(length(which(x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(geneName)),]$Variant_Classification %in% c("Frame_Shift_DEL", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation")))>0, 1, 0)))
  likely <- ifelse(polyphen == 1 | sift == 1, 1, 0)
  return(list(likely, definite))
}


#Mut - select one
mut_value <- function(muts, gene, selection, amino){
  if(selection == "Select my own mutations"){
    gene_mut <- unlist(lapply(muts, function(x) {
      ifelse(length(which(paste0(x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$HUGO_SYMBOL, ": ",
        gsub("p.", "", x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$HGVSp_Short)) %in% amino))>0, 1, 0)
    }))
  } else if (selection == "Inactivating"){
    gene_mut <- calculateMutType(muts, gene)[[2]]
  } else if (selection == "Likely inactivating"){
    gene_mut <- calculateMutType(muts, gene)[[1]]
  } else if (selection == "None"){
    gene_mut<- NA
  }
  return(gene_mut)
}

#CNA - select one
cna_value <- function(cna, gene, selection){
 if(length(selection) > 0){
    gene_cna <- unlist(lapply(cna, function(x) {
      ifelse(length(which(paste0(x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$HUGO_SYMBOL, ": ",
                                 x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$ALTERATION) %in% selection))>0, 1, 0)
    })) 
    } else {
      gene_cna <- NA
    }
  return(gene_cna)
}

#Fusion - select one
fusion_value <- function(fusion, gene, selection){
  if(length(selection)>0){
    gene_fusion <- unlist(lapply(fusion, function(x) {
      ifelse(length(which(paste0(x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$HUGO_SYMBOL, ": ",
                                 x[which(as.character(unlist(x$HUGO_SYMBOL)) %in% unlist(gene)),]$MUTATION_EFFECT) %in% selection))>0, 1, 0)
    }))
    } else {
    gene_fusion <- NA
  }
  return(gene_fusion)
}

