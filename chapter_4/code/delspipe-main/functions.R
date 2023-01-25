# Filter enzymatic deletions. 
filtEnzDels <- function(dels){
        enzDels <- dels[!is.na(dels$EC_number), ]
        return(enzDels)
}

# Binarize deletion files. 
binarizeDels <- function(DM, thrshld){
        nums <- DM[, sapply(DM, class) == "numeric"]
        gInf <- DM[, sapply(DM, class) != "numeric"]
        bin <- nums
        bin[bin < thrshld] <- 0
        bin[bin >= thrshld] <- 1
        binOut <- cbind.data.frame(gInf, bin)
        return(binOut)
}