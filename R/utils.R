## just simple utility function to test for empty/not set column names
noColNames <- function(x){
    CN <- colnames(x)
    if(all(CN %in% paste("V", 1:ncol(x), sep=""))){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

## takes a pedigree (data.frame, pedigree or pedigreeList) as input and generates a
## valid pedigree for the submitted ids (if no ids are specified it returns the ped).
## basically, the function subsets the full pedigree to the individuals ids and
## recursively adds parents if needed.
buildPedigree <- function(ped, ids){
    ## check pedigree and eventually add missing parents...
    if(missing(ids))
        return(ped)
    ## want to have a data.frame or matrix...
    if(is(ped, "pedigree") | is(ped, "pedigreeList"))
        ped <- ped2df(ped)
    if(!all(ids %in% rownames(ped)))
        stop("Some of the submitted ids are not in the pedigree!")
    subfam <- ped[as.character(ids), ]
    ## now recursively add parents, if needed...
    parents <- as.character(unique(c(subfam$father, subfam$mother)))
    ## parents[is.na(parents)] <- 0
    ## parents <- parents[parents!=0]
    parents <- parents[!is.na(parents)]
    while(!all(parents %in% as.character(subfam$id))){
        subfam <- ped[unique(c(as.character(subfam$id), parents)), ]
        parents <- as.character(unique(c(subfam$father, subfam$mother)))
        ## parents[is.na(parents)] <- 0
        ## parents <- parents[parents!=0]
        parents <- parents[!is.na(parents)]
    }
    return(subfam)
}


doGetSiblings <- function(ped, id=NULL){
    if(is.null(id))
        stop("At least one if has to be specified!")
    parents <- doGetAncestors(ped, id=id, maxlevel=1)
    siblings <- doGetChildren(ped, id=parents, maxlevel=1)
    return(siblings)
}

## gets all the ancestors for (one or more) id(s)
doGetAncestors <- function(ped, id=NULL, maxlevel=3){
    if(is.null(id))
        stop("At least one id has to be specified!")
    allids <- NULL
    for(i in 1:maxlevel){
        ## get the parents of all the allids
        newids <- unlist(
            ped[as.character(ped$id) %in% id, c("father", "mother")],
            use.names=FALSE
        )
        zeros <- is.na(newids)
        ##zeros <- newids == founder
        if(all(zeros) | length(newids)==0)
            break
        newids <- newids[!zeros]
        id <- newids
        allids <- c(allids, newids)
    }
    return(as.character(unique(allids)))
}

## gets all the children for (one or more) id(s)
## note: these are only children in direct blood line.
doGetChildren <- function(ped, id=NULL, maxlevel=16){
    if(is.null(id))
        stop("At least one id has to be specified!")
    allids <- NULL
    for(i in 1:maxlevel){
        ## get all the children for the ids
        newids <- ped[ped$father %in% id | ped$mother %in% id, "id"]
        if(length(newids)==0)
            break
        id <- newids
        allids <- c(allids, newids)
    }
    return(as.character(unique(allids)))
}

## check for the submitted id(s) whether a mate (spouse) is
## available (in column mother or father of the pedigree ped)
## and if yes add them to the returned list of ids.
## returned value contains
## eventually missing mates.
doGetMissingMate <- function(ped, id){
    if(is.null(id))
        stop("Al least one id has to be specified")
    newids <- c(ped[ped$father %in% id, "mother"],
                ped[ped$mother %in% id, "father"])
    return(unique(newids))
}


writePed <- function(ped, file){
    Sexes <- as.numeric(sanitizeSex(ped$sex))
    Sexes[is.na(Sexes)] <- 0
    ped$sex <- Sexes
    write.table(ped, file=file, quote=FALSE, sep="\t", row.names=FALSE,
                col.names=FALSE)
}

## creates a data.frame (or matrix) for a pedigree...
ped2df <- function(ped){
    ped <- do.call(cbind, ped)
    if(!any(colnames(ped)=="famid"))
        ped <- cbind(ped, famid=rep(1, nrow(ped)))
    ped <- ped[ , c("famid", "id", "findex", "mindex", "sex")]
    ## FIXME!!! Can I get NAs as findex too?
    ## have to rematch the fater id and mother id...
    notZ <- which(ped[, "findex"] > 0)
    ped[notZ, "findex"] <- ped[ped[notZ, "findex"], "id"]
    notZ <- which(ped[, "findex"] > 0)
    ped[notZ, "mindex"] <- ped[ped[notZ, "mindex"], "id"]
    colnames(ped) <- .PEDCN
    ped <- data.frame(ped)
    rownames(ped) <- ped$id
    return(ped)
}


df2ped <- function(ped){
    ped <- checkPedCol(ped)
    if(!any(colnames(ped)=="family"))
        ped <- cbind(ped, family=rep(1, nrow(ped)))
    ## if there is only one family, remove that one ;)
    if(length(unique(ped$family)) == 1)
        ped <- ped[, colnames(ped)!="family"]
    ##
    if(any(colnames(ped) == "affected")){
        if(any(colnames(ped) == "family")){
            peddi <- kinship2::pedigree(id=ped$id, dadid=ped$father, momid=ped$mother,
                                        famid=ped$family, sex=as.numeric(ped$sex),
                                        affected=ped$affected)
        }else{
            peddi <- kinship2::pedigree(id=ped$id, dadid=ped$father, momid=ped$mother,
                                        sex=as.numeric(ped$sex),
                                        affected=ped$affected)
        }
    }else{
        if(any(colnames(ped) == "family")){
            peddi <- kinship2::pedigree(id=ped$id, dadid=ped$father, momid=ped$mother,
                                        famid=ped$family, sex=as.numeric(ped$sex))
        }else{
            peddi <- kinship2::pedigree(id=ped$id, dadid=ped$father, momid=ped$mother,
                                        sex=as.numeric(ped$sex))
        }
    }
    return(peddi)
}


## finds related individuals, i.e. individuals sharing the same blood line.
## should we call this "shareBloodLine" or "shareKinship"
## ped has to be the data.frame that we get with the pedigree method.
doShareKinship <- function(ped, kin, id){
    if(missing(ped) & missing(kin))
        stop("ped or kin has to be submitted!")
    if(missing(id))
        stop("At least one id has to be submitted!")
    id <- as.character(id)
    ## calculate kin from ped
    if(missing(kin)){
        kin <- kinship(id=ped$id, momid=ped$mother, dadid=ped$father,
                    sex=ped$sex)
    }
    if(!all(id %in% colnames(kin)))
        stop("id has to be present in ped or kin!")
    ## subset to id (which can be more than one)
    kin <- kin[id, , drop=FALSE]
    related <- colnames(kin)[colSums(kin) > 0]
    return(related)
}

####============================================================
##  doGetFounders
##
##  Returns the ids of the founders in the specified pedigree.
####------------------------------------------------------------
doGetFounders <- function(ped){
    founders <- ped[is.na(ped$mother) & is.na(ped$father), "id"]
    return(founders)
}

####============================================================
##  doGetSingletons
##
##  Return the ids of those founders that have no child in the
##  pedigree.
####------------------------------------------------------------
doGetSingletons <- function(ped){
    founders <- ped[is.na(ped$mother) & is.na(ped$father), , drop=FALSE]
    if(nrow(founders) > 0){
        childl <- !((founders$id %in% ped$father) | (founders$id %in% ped$mother))
        if(any(childl)){
            return(founders[childl, "id"])
        }
    }
    return(character())
}


## gets all parents.
getParents <- function(ped, id){
    if(is(ped, "pedigree") | is(ped, "pedigreeList"))
        ped <- ped2df(ped)
    ped <- checkPedCol(ped)
    if(missing(id))
        stop("id is missing!")
    if(is.null(rownames(ped)))
        rownames(ped) <- ped[, "id"]
    id <- as.character(id)
    if(!is.data.frame(ped))
        ped <- as.data.frame(ped)
    return(ped[id, c("father", "mother"), drop=FALSE])
}


## just a simple helper function that checks that the ped parameter
## has the required colnames, is a data.frame and has rownames set.
checkPedCol <- function(ped){
    if(missing(ped))
        stop("ped is required!")
    if(!all(c("id", "father", "mother") %in% colnames(ped)))
        stop("id, father, mother are required columns in ped.")
    if(!is.data.frame(ped))
        ped <- as.data.frame(ped)
    ## add a family column if not present...
    if(!any(colnames(ped) == "family"))
        ped <- cbind(ped, family=rep(1, nrow(ped)))
    rownames(ped) <- as.character(ped$id)
    return(ped)
}

## Sanitize the pedigree the way we want: sex should be sanitized (if present),
## father and mother column should contain NA for founders instead of "" or 0.
sanitizePed <- function(ped){
    if(any(colnames(ped) == "sex"))
        ped$sex <- sanitizeSex(ped$sex)
    Father <- ped$father
    if(is.factor(Father))
        Father <- as.character(Father)
    if(is.numeric(Father))
        Father[which(Father == 0)] <- NA
    if(is.character(Father))
        Father[which(Father == "")] <- NA
    Mother <- ped$mother
    if(is.factor(Mother))
        Mother <- as.character(Mother)
    if(is.numeric(Mother))
        Mother[which(Mother == 0)] <- NA
    if(is.character(Mother))
        Mother[which(Mother == "")] <- NA
    ped$father <- Father
    ped$mother <- Mother
    return(ped)
}

## find the subPedigree (smallest pedigree) includin the individuals specified
## with id and all eventually needed additionals
subPedigree <- function(ped, id=NULL, all=TRUE){
    CL <- class(ped)
    ## handling for pedigree and pedigreeList classes: unfortunately these are not
    ## correct classes, otherwise I could define methods instead of a function...
    if(is(ped, "pedigree") | is(ped, "pedigreeList")){
        ped <- ped2df(ped)
    }
    ped <- checkPedCol(ped)
    if(is.null(id))
        stop("id is missing!")
    id <- as.character(id)
    subgr <- connectedSubgraph(ped2graph(ped), nodes=id, mode="all", all.nodes=all,
                               ifnotfound=NA)
    if(!is(subgr, "igraph")){
        stop("No common pedigree for the submitted individuals found!")
    }
    inds <- names(V(subgr))
    ## get missing mates:
    inds <- unique(c(inds, doGetMissingMate(ped, id=inds)))
    subped <- ped[inds, ]
    subped <- removeSingletons(subped)
    ## fix founders.
    subped[!(subped$father %in% subped$id), "father"] <- NA
    subped[!(subped$mother %in% subped$id), "mother"] <- NA
    if(CL == "pedigree" | CL == "pedigreeList")
        subped <- df2ped(subped)
    return(subped)
}




## find the (eventually smallest) connected subgraph of all specified
## nodes.
## graph: igraph object.
## nodes: node (vertex) names.
## mode: to calculate the smallest distance.
## all.nodes: if all nodes should be present in the resulting graph.
## ifnotfound: by default, the function throws an error if no connected subgraph is
##             found, or, if all.nodes=TRUE, not all input nodes are in the connected
##             subgraph. If ifnotfound is specified, its value is returned instead of
##             an error.
connectedSubgraph <- function(graph, nodes, mode="all", all.nodes=TRUE, ifnotfound){
    if(!is(graph, "igraph"))
        stop("graph has to be an igraph object!")
    mode <- match.arg(mode, c("all", "in", "out"))
    if(missing(nodes) | length(nodes) < 2)
        stop("At least two nodes have to be submitted!")
    nodes <- as.character(nodes)
    ## check if nodes are present in the graph
    nono <- !(nodes %in% names(V(graph)))
    if(any(nono)){
        missingNodes <- nodes[nono]
        if(all.nodes)
            stop("With argument all.nodes=TRUE all nodes have to be present in the graph! Nodes ",
                 paste(missingNodes, collapse=", "), " not present in the original graph.")
        nodes <- nodes[!nono]
        warning(paste0("Nodes ", paste(missingNodes, collapse=", "),
                       " removed as they are not present in graph! "))
        if(length(nodes) < 2)
            stop("At least two nodes have to be present in the graph!")
    }
    ## OK, now we can start...
    ## 1) iterate over all nodes and calculate the shortest path between all.
    ## 2) get all nodes along the shortest path and add to the needed nodes.
    needNodes <- character()
    ## two versions: if mode="all" the directionality does not matter, thus we
    ## can run a faster version. for the other modes really calculate distances
    ## between all (both directions).
    if(mode == "all"){
        for(i in 1:(length(nodes)-1)){
            suppressWarnings(
                paths <- shortest_paths(graph, from=nodes[i],
                                        to=nodes[(i+1):length(nodes)],
                                        mode=mode)
            )
            needNodes <- unique(c(needNodes, unlist(lapply(paths$vpath, names))))
        }
    }else{
        for(i in 1:length(nodes)){
            suppressWarnings(
                paths <- shortest_paths(graph, from=nodes[i],
                                        to=nodes[-i],
                                        mode=mode)
            )
            needNodes <- unique(c(needNodes, unlist(lapply(paths$vpath, names))))
        }
    }
    if(length(needNodes)==0){
        if(missing(ifnotfound)){
            stop("Could not find any connected subgraph!")
        }
        warning("Could not find any connected subgraph!")
        return(ifnotfound)
    }
    ## check the subgraph.
    if(!all(nodes %in% needNodes)){
        if(all.nodes){
            if(missing(ifnotfound)){
                stop("Could not find any connected subgraph!")
            }
            warning("Could not find any connected subgraph!")
            return(ifnotfound)
        }else{
            warning(paste0("Nodes ",
                           paste(nodes[!(nodes %in% needNodes)]),
                           " not present in the returned graph!")
                    )
        }
    }
    ## finally build the graph...
    ## subgraph is a little tedious as it requires submitting of vertex INDICES!
    return(induced_subgraph(graph, vids=needNodes))
}

## search for a common ancestor couple in a pedigree.
## ped: data.frame or pedigree or pedigreeList.
doGetCommonAncestor <- function(ped, id, method="min.dist"){
    method <- match.arg(method, c("min.dist", "smallest.mean.dist"))
    if(missing(id))
        stop("id has to be specified")
    if(length(id) <= 1)
        stop("The minimum number of ids has to be 2.")
    if(is(ped, "pedigree") | is(ped, "pedigreeList"))
        ped <- ped2df(ped)
    ped <- checkPedCol(ped)
    ## check if id are in ped:
    id <- as.character(id)
    if(!all(id %in% rownames(ped))){
        stop("all of the sumitted ids should be present in the pedigree!")
    }
    ## check if the ids are in the same family
    if(any(colnames(ped)=="family")){
        if(length(unique(ped[id, "family"])) > 1){
            warning("The provided ids are from different families; still",
                    " trying to find a common ancestor but I'll probably fail.")
        }else{
            ## subset the pedigree to this family...
            ped <- ped[ped$family == ped[id, "family"][1], , drop=FALSE]
        }
    }
    ## OK; done with input parameter checking.
    ## transform the pedigree into a graph
    igr <- ped2graph(ped)
    ## "in": go from children -> parents
    ## LLLL TODO can do faster!!!
    din <- distances(igr, mode="in")
    ## use the ROWS in that.
    subD <- din[id, , drop=FALSE]
    subD <- subD[, !apply(subD, MARGIN=2, function(x){
        return(all(infOrZero(x)))
    }) , drop=FALSE]
    ## check if we've got one common ancestor:
    if(all(is.infinite(colSums(subD)))){
        warning("Did not find any common ancestor!")
        return(NA)
    }
    ## remains to find the "one and only" (well, usually we expect a couple).
    subD <- subD[, !is.infinite(colSums(subD)), drop=FALSE]
    ## now get the ones that are closest:
    if(method == "min.dist"){
        dists <- apply(subD, MARGIN=2, min)
        ancs <- names(dists)[dists==min(dists)]
    }
    if(method == "smallest.mean.dist"){
        dists <- colMeans(subD)
        ancs <- names(dists)[dists==min(dists)]
    }
    ## at last trying to get the mate... if any; usually, ancs should
    ## already correspond to the ancestor couple.
    ancs <- c(ancs, doGetMissingMate(ped, ancs))
    return(unique(ancs))
}

## returns true if a value is 0 or infinite
infOrZero <- function(x){
    return((is.infinite(x) | x == 0))
}

## transform a pedigree into a (directed) igraph
## direction is always parent->child.
## this time each individual is treated as a single node!
ped2graph <- function(ped){
    if(is(ped, "pedigree") | is(ped, "pedigreeList"))
        ped <- ped2df(ped)
    ped <- checkPedCol(ped)
    ## define end/starting points
    ped[!ped[, "father"] %in% ped$id, "father"] <- NA
    ped[!ped[, "mother"] %in% ped$id, "mother"] <- NA
    eds <- apply(ped[, c("id", "father", "mother")], MARGIN=1, function(x){
        retval <- c()
        if(!is.na(x["father"]))
            retval <- x[c("father", "id")]
        if(!is.na(x["mother"]))
            retval <- c(retval, x[c("mother", "id")])
        return(retval)
    })
    ## what the heck... I need character vertex names!
    Edges <- unlist(eds, use.names=FALSE)
    igr <- graph(edges=as.character(Edges))
    return(igr)
}


## takes a pedigree as input and tries to identify the founders
## in case there are more than two founder pairs (i.e. founders with the same, largest
## number of offspring generations), it returns only one of them.
doFindFounders <- function(ped, family){
    if(!any(colnames(ped) == "family"))
        ped <- data.frame(ped, family=rep(1, nrow(ped)), stringsAsFactors=FALSE)
    families <- unique(as.character(ped$family))
    if(missing(family))
        family <- NULL
    if(is.null(family)){
        family <- families[1]
        warning(paste0("No family specified; using the first one in the pedigree: ",
                       family, "."))
    }
    ped <- ped[ped$family == family, , drop=FALSE]
    ## OK, now we have the pedigree of a single family... try to find founders.
    founders <- ped$id[which(is.na(ped$father) & is.na(ped$mother))]
    ## determine for each founder the number of generations of its children, grandchildren etc.
    founderGen <- doCountGenerations(ped, id=founders, direction="down")
    founders <- names(founderGen)[founderGen == max(founderGen)][1]
    founders <- c(founders, doGetMissingMate(ped, founders))
    return(unique(founders))
}


## calculate for each id the number of generations up- or down.
doCountGenerations <- function(ped, id=NULL, direction="down"){
    if(is.null(id))
        stop("At least one id has to be specified!")
    direction <- match.arg(direction, c("up", "down"))
    if(direction == "up")
        walkFun <- doGetAncestors
    if(direction == "down")
        walkFun <- doGetChildren
    res <- sapply(id, function(x){
        generation <- 0
        ids <- x
        repeat{
            ids <- walkFun(ped, id=ids, maxlevel=1)
            if(length(ids) == 0){
                return(generation)
            }else{
                generation <- generation + 1
            }
        }
    })
    names(res) <- id
    return(res)
}


## returns a list of generations for each family in the pedigree.
## generations are relative to the founders.
doEstimateGenerationsFor2 <- function(ped, family=NULL){
    if(!is.null(family)){
        ped <- ped[ped$family %in% family, ]
        if(nrow(ped) == 0)
            stop("The submitted family id is not present in the pedigree!")
    }
    ## OK, now loop through the families...
    pedList <- split(ped, ped$family)
    Gens <- lapply(pedList, function(currentPed){
        founders <- doFindFounders(currentPed, family=currentPed[1, "family"])
        generation <- doGetGenerationFrom2(currentPed, founders[1])
        return(generation)
    })
    return(Gens)
}


getAllSibsNMates <- function(ped, id){
    ## recursively fetch siblings and mates
    sibsNmates <- id
    repeat{
        sibs <- unique(c(id, doGetSiblings(ped, id=id)))
        sibsNmates <- unique(c(doGetMissingMate(ped, id=sibs), sibs))
        if(length(sibsNmates) != length(id)){
            id <- sibsNmates
        }else{
            break
        }
    }
    return(sibsNmates)
}

doGetGenerationFrom2 <- function(ped, id, generation=0){
    if(missing(id))
        stop("id has to be specified!")
    if(missing(ped))
        stop("ped has to be specified!")
    id <- as.character(id)
    if(length(id) > 1){
        id <- id[1]
        warning("Length of 'id' is larger 1. Using only the first id in 'id': ", id)
    }
    if(!any(ped$id %in% id))
        stop("Can not find id in the pedigree!")
    famid <- ped[ped$id == id, "family"]
    ped <- ped[ped$family == famid, , drop=FALSE]
    ## first define the connected!!! pedigree.
    pedg <- ped2graph(ped)
    if(clusters(pedg)$no > 1){
        Cl <- clusters(pedg)
        warning("Found ", Cl$no, " connected sub-graphs in the pedigree for family ", famid,
                "! Will restrict the analysis to the subgraph containing individual ", id, ".")
        vertices <- names(Cl$membership)[Cl$membership == Cl$membership[id]]
        pedg <- induced_subgraph(pedg, vids=vertices)
    }
    individuals <- names(V(pedg))
    gens4individuals <- rep(NA, length(individuals))
    names(gens4individuals) <- individuals
    if(!any(names(gens4individuals) == id))
        stop(paste0("Individual ", id, " is not part of the connected pedigree!"))
    gens4individuals[id] <- 0
    fromId <- id
    Tested <- rep(FALSE, length(gens4individuals))    ## that helps to prevent testing the same individual several times...
    repeat{
        ## get the generation of the individuals from which we start
        fromIdBool <- individuals %in% fromId
        Tested[fromIdBool] <- TRUE
        currentGens <- gens4individuals[fromIdBool]
        currentGens <- split(currentGens, currentGens)
        Res <- lapply(currentGens, allids=individuals, FUN=function(z, allids){
            ## get the ancestors and set the generation for them
            ancs <- getAncestors(ped, id=names(z), max.generations=1)
            if(length(ancs) > 0){
                ancret <- cbind(which(allids %in% ancs), (z[1] - 1))
            }else{
                ancret <- NULL
            }
            ## get the children
            kids <- getChildren(ped, id=names(z), max.generation=1)
            if(length(kids) > 0){
                kidret <- cbind(which(allids %in% kids), (z[1] + 1))
            }else{
                kidret <- NULL
            }
            return(rbind(ancret, kidret))
        })
        Res <- do.call(rbind, Res)
        if(nrow(Res) == 0){
            warning("What the heck... stopped before we got all?")
            break
        }
        ## don't want to search the same id several times... if Tested is TRUE, don't use them again.
        toTestIdx <- Res[, 1][!Tested[Res[, 1]]]
        fromId <- individuals[toTestIdx]
        gens4individuals[Res[, 1]] <- Res[, 2]
        if(sum(is.na(gens4individuals)) == 0)
            break
    }
    ## make a larger generation vector.
    allgens <- rep(NA, nrow(ped))
    names(allgens) <- ped$id
    allgens[names(gens4individuals)] <- gens4individuals
    return(allgens)
}

##********************************************************************
##
##   Utility functions
##
##********************************************************************
sanitizeSex <- function(x){
    if(is.numeric(x)){
        ## require values 1, 2, NA
        numrange <- unique(x)
        if(!all(c(1, 2) %in% numrange))
            stop(paste0("If column sex is numeric it has to contain 1 (=male) ",
                        "and 2 (=female) values! Unknown or missing can be encoded",
                        " by NA or any other number than 1 or 2."))
        news <- rep(NA, length(x))
        news[which(x == 1)] <- "M"
        news[which(x == 2)] <- "F"
        return(factor(news, levels=c("M", "F")))
    }
    if(is.factor(x))
        x <- as.character(x)
    if(is.character(x)){
        ## remove all white spaces:
        x <- gsub(x, pattern="\\s", replacement="")
        ## that's tricky...kind of...
        ## just get the first character; want to have M and F...
        x <- substr(toupper(x), 1, 1)
        chars <- unique(x)
        if(length(chars) == 2){
            if(!all(c("M", "F") %in% chars))
                stop(paste0("Male and female individuals have to be represented",
                            " by 'M' and 'F', respectively!"))
        }
        if(length(chars) == 1){
            if(!any(c("M", "F") %in% chars))
                stop(paste0("Male and female individuals have to be represented",
                            " by 'M' and 'F', respectively!"))
        }
        if(length(chars[!is.na(chars)]) > 2)
            warning(paste0("All characters in column sex other than M and F",
                           " (or male, female etc) are set to NA!"))
        news <- rep(NA, length(x))
        news[which(x == "M")] <- "M"
        news[which(x == "F")] <- "F"
        return(factor(news, levels=c("M", "F")))
    }
}


#####
#### NOTE:
##   isn't that the same as removeSingletons???
##
####
## remove non-connected individuals. We're assuming that ped is the pedigree
## of a SINGLE family, as we expect to get ONE graph. in case we've got 2, we're
## using the largest connected one.
## solveMultiGraph: what should we do if we get multiple graphs?
##    + use.all: just use all, i.e. return vertix names of all.
##    + use.largest: just return vertices of the largest connected sub-graph.
doPrunePed <- function(ped, addMissingMates=FALSE, solveMultiGraph="use.all"){
    .Deprecated("removeSingletons")
    solveMultiGraph <- match.arg(solveMultiGraph, c("use.all", "use.largest"))
    pedg <- ped2graph(ped)
    if(solveMultiGraph == "use.largest"){
        if(clusters(pedg)$no > 1){
            Cl <- clusters(pedg)
            warning("Found ", Cl$no,
                    " connected sub-graphs in the pedigree!",
                    " Will restrict the analysis to the largest connected pedigree",
                    " in that family.")
            vertices <- names(Cl$membership)[Cl$membership == as.numeric(sort(table(Cl$membership),
                                                                              decreasing=TRUE))[1]]
            pedg <- induced_subgraph(pedg, vids=vertices)
        }
    }
    keepIds <- names(V(pedg))
    if(addMissingMates){
        missMates <- doGetMissingMate(ped, id=keepIds)
        keepIds <- unique(c(missMates, keepIds))
    }
    notThere <- !(as.character(ped$id) %in% keepIds)
    ped <- ped[!notThere, , drop=FALSE]
    if(any(notThere))
        warning("Had to remove ", sum(notThere),
                " un-connected individuals from the pedigree.")
    return(ped)
}


## a singleton a.k.a childless founder is an individual which id is not present in the
## father or mother column of the pedigree and which has 0 in mother and father
.removeSingletons <- function(ped){
    message("Removing singletons...", appendLF=FALSE)
    ped <- checkPedCol(ped)
    ped <- sanitizePed(ped)
    childlF <- doGetSingletons(ped)
    if(length(childlF) > 0){
        ped <- ped[!(ped$id %in% childlF), ]
        message(" ", length(childlF), " removed.")
    }else{
        message(" none present.")
    }
    return(ped)
}

