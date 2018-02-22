## pedFile <- system.file("txt/Large.ped", package="FamAgg")
## cliqFile <- system.file("txt/Large-cliques.txt", package="FamAgg")
## ped <- read.table(pedFile, sep=";")
## cliqDf <- read.table(cliqFile, sep="\t")
## cliqs <- cliqDf[, 1]
## names(cliqs) <- cliqDf[, 2]
## ## building the pedigree
## fad <- FAData(ped[, 1:5])
## affs <- ped[, 6]
## affs[affs==2] <- 1
## names(affs) <- ped[, 2]


## ## OK, 0.0.6
## ## new:
## ## compare families with cliques from minnbreast_default_cliques.txt
## test_probability <- function(){
##     ## pedFile <- system.file("txt/Large.ped", package="FamAgg")
##     ## cliqFile <- system.file("txt/Large-cliques.txt", package="FamAgg")
##     ## ped <- read.table(pedFile, sep=";")
##     ## cliqDf <- read.table(cliqFile, sep="\t")
##     ## cliqs <- cliqDf[, 1]
##     ## names(cliqs) <- cliqDf[, 2]
##     ## ## building the pedigree
##     ## fad <- FAData(ped[, 1:5])
##     ## affs <- ped[, 6]
##     ## affs[affs==2] <- 1
##     ## names(affs) <- ped[, 2]
##     ## performing the test
##     if (.Platform$OS.type == "unix") {
##         far <- probabilityTest(fad, trait=affs,
##                                cliques=cliqs, nsim=500)
        
##         ## check if the cliques is the same as cliqs
##         farCliqs <- cliques(far, na.rm=TRUE)
##         checkEquals(as.character(farCliqs[names(cliqs)]), as.character(cliqs))
        
##         ## shuffle cliqs and replace; check if we get the correct one back.
##         cliques(far) <- sample(cliqs, size=length(cliqs))
##         farCliqs <- cliques(far, na.rm=TRUE)
##         checkEquals(as.character(farCliqs[names(cliqs)]), as.character(cliqs))
##         far <- runSimulation(far, nsim=500)
        
##         head(result(far))
        
##         ## testing to subset the object... which is not supported
##         checkException(far[1:10, ])
##     }
##     return(TRUE)
## }


## test_plot_probability <- function(){
##     do.plot <- FALSE
##     ## pedFile <- system.file("txt/Large.ped", package="FamAgg")
##     ## cliqFile <- system.file("txt/Large-cliques.txt", package="FamAgg")
##     ## ped <- read.table(pedFile, sep=";")
##     ## cliqDf <- read.table(cliqFile, sep="\t")
##     ## cliqs <- cliqDf[, 1]
##     ## names(cliqs) <- cliqDf[, 2]
##     ## ## building the pedigree
##     ## fad <- FAData(ped[, 1:5])
##     ## affs <- ped[, 6]
##     ## affs[affs==2] <- 1
##     ## names(affs) <- ped[, 2]
##     ## performing the test
##     if (.Platform$OS.type == "unix") {
##         far <- probabilityTest(fad, trait=affs,
##                                cliques=cliqs, nsim=500)
##         res <- result(far)
##         id <- res[1, "group_id"]
##         ped <- buildPed(far, id=id, prune=TRUE)
##         allCliqs <- cliques(far)
##         clique <- names(allCliqs)[which(allCliqs == as.character(id))]
##         ## plot the pedigree highlighting the clique.
##         ## plotPed(far, id=id, filename="~/Desktop/test.pdf", prune=TRUE)
##         switchPlotfun("ks2paint")
##         plotPed(far, id=id, prune=TRUE, device="plot")
##         ## the large huge busy pedigree.
##         if(do.plot){
##             plotPed(far, id=id, prune=FALSE, device="plot")
##         }
##         id <- res[2, "group_id"]
##         ped <- buildPed(far, id=id, prune=TRUE)
##         clique <- names(allCliqs)[which(allCliqs == as.character(id))]
##         all(clique %in% ped$id)
##         plotPed(far, id=id, prune=TRUE, device="plot")
##         ##switchPlotfun()
##     }
##     return(TRUE)
## }

## test_buildped_probability <- function(){
##     ## pedFile <- system.file("txt/Large.ped", package="FamAgg")
##     ## cliqFile <- system.file("txt/Large-cliques.txt", package="FamAgg")
##     ## ped <- read.table(pedFile, sep=";")
##     ## cliqDf <- read.table(cliqFile, sep="\t")
##     ## cliqs <- cliqDf[, 1]
##     ## names(cliqs) <- cliqDf[, 2]
##     ## ## building the pedigree
##     ## fad <- FAData(ped[, 1:5])
##     ## affs <- ped[, 6]
##     ## affs[affs==2] <- 1
##     ## names(affs) <- ped[, 2]
##     ## performing the test
##     if (.Platform$OS.type == "unix") {
##         far <- probabilityTest(fad, trait=affs,
##                                cliques=cliqs, nsim=500)
##         res <- result(far)
##         id <- res[1, "group_id"]
##         ## get the clique
##         clique <- names(cliqs)[which(as.character(cliqs) == as.character(id))]
##         ## build the full pedigree
##         FullPed <- buildPed(far, id=id)
##         ## check if we've got all of the clique in the pedigree...
##         checkEquals(length(clique), sum(clique %in% FullPed$id))
##         ## make a small subset pedigree:
##         SmallPed <- buildPed(far, id=id, prune=TRUE)
##         ## check if we've got all of the clique in the pedigree...
##         checkEquals(length(clique), sum(clique %in% SmallPed$id))
##     }
##     return(TRUE)
## }





