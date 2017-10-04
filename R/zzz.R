### setting options and stuff...
.setOptions <- function(){
    FA <- list()
    FA$haplopainter <- system.file("perl/HaploPainter1.043.pl",
                                   package="FamAgg")
    options("FamAgg"=FA)
}


.checkHaplopaintRequirements <- function(messagefun=message){
    pedfile <- system.file("txt/minnbreastsub.txt", package="FamAgg")
    mbsub <- read.table(pedfile, sep="\t", header=TRUE)
    PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
    colnames(PedDf) <- .PEDCN
    PedDf$sex <- sanitizeSex(PedDf$sex)
    suppressMessages(
        fad <- FAData(PedDf)
    )
    ped <- buildPed(fad, id=6)
    pedDf <- buildHaplopaintDataframe(individual=ped$id, father=ped$father,
                                      mother=ped$mother, gender=ped$sex)
    tmpf <- tempfile()
    outfile <- paste0(tmpf, ".pdf")
    write.table(pedDf, tmpf, row.names=FALSE, col.names=FALSE,
                sep="\t", quote=FALSE)
    plotcall <- paste0("perl ", options()$FamAgg$haplopaint, " -b -pedfile ",
                       tmpf, " -pedformat csv -outfile ", outfile,
                       " -family 1 -bgcolor \\#ffffff -outformat pdf")
    res <- tryCatch(system(plotcall), error=function(e){
    })
    if(res == 0){
        messagefun("OK")
        return("haplopaint")
    }else{
        messagefun(paste0("FAIL\nHaplopainter not working! For requirements ",
                          "and installation see the package vignette.\nWill ",
                          "use internal plotting functions instead."))
        return("ks2paint")
    }
}

.onLoad <- function(libname, pkgname){
    .setOptions()
}

.onAttach <- function(libname, pkgname){
    .setOptions()
    FA <- getOption("FamAgg")
    FA$plotfun <- "ks2paint"
    options("FamAgg"=FA)
}

## maximal supported pedigree size for the gap package:
GAP_MAX_CLIQUE_SIZE <- 22

