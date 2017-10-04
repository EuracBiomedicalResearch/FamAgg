## here we define the pedigree plotting functions.
switchPlotfun <- function(method){
    if(!missing(method)){
        method <- match.arg(method, c("ks2paint", "haplopaint"))
    }else{
        if(options()$FamAgg$plotfun=="ks2paint"){
            method <- "haplopaint"
        }else if(options()$FamAgg$plotfun=="haplopaint"){
            method <- "ks2paint"
        }
    }
    if(method=="haplopaint"){
        ## check if that is possible at all.
        method <- .checkHaplopaintRequirements()
    }
    FA <- getOption("FamAgg")
    FA$plotfun <- method
    message(paste0("Plotting backend: ", method, "."))
    options("FamAgg"=FA)
}

##########################
## the main plotting function that dispatches on options()$FamAgg$plotfun
## this function takes a lot of input...
doPlotPed <- function(family=NULL, individual=NULL, father=NULL, mother=NULL,
                      gender=NULL, affected=NULL, is.deceased=NULL,
                      is.sab.or.top=NULL, is.proband=NULL, is.adopted=NULL,
                      are.twins=NULL, are.consanguineous=NULL,
                      text.inside.symbol=NULL, text.beside.symbol=NULL,
                      text1.below.symbol=NULL, text2.below.symbol=NULL,
                      text3.below.symbol=NULL, text4.below.symbol=NULL,
                      filename=NULL, main=NULL, device="plot", res=600, ...){
    plotfun <- options()$FamAgg$plotfun
    plotfun <- match.arg(plotfun, c("ks2paint", "haplopaint"))
    if(plotfun == "haplopaint"){
        if(device == "plot"){
            device="pdf"
            warning("haplopaint does not support device='plot', changing ",
                    "to 'pdf'.")
        }
    }
    ## build the argument list.
    argList <- list(family=family, individual=individual, father=father,
                    mother=mother, gender=gender, affected=affected,
                    is.deceased=is.deceased, is.sab.or.top=is.sab.or.top,
                    is.proband=is.proband, is.adopted=is.adopted,
                    are.twins=are.twins,
                    are.consanguineous=are.consanguineous,
                    text.inside.symbol=text.inside.symbol,
                    text.beside.symbol=text.beside.symbol,
                    text1.below.symbol=text1.below.symbol,
                    text2.below.symbol=text2.below.symbol,
                    text3.below.symbol=text3.below.symbol,
                    text4.below.symbol=text4.below.symbol,
                    filename=filename, main=main, device=device,
                    res=res)
    addArgs <- list(...)
    ## kick out arguments we know might be passed, but are not needed
    ## further down...
    if(length(addArgs) > 0){
        addArgs <- addArgs[!(names(addArgs) %in%
                             c("prune", "max.generations.up",
                               "max.generations.down"))]
    }
    argList <- c(argList, addArgs)
    ## call the function.
    invisible(do.call(plotfun, args=argList))
}

##########################
## ks2paint uses kinship2 to plot pedigrees.
## text2.below.symbol will be plotted on the top left corner.
## text3.below.symbol will be plotted on the top right corner.
## text1.below.symbol is expected to be the age... add that below the ID.
## is.proband: highlight ids in red.
ks2paint <- function(family=NULL, individual=NULL, father=NULL, mother=NULL,
                     gender=NULL, affected=NULL, is.deceased=NULL,
                     is.sab.or.top=NULL, is.proband=NULL, is.adopted=NULL,
                     are.twins=NULL, are.consanguineous=NULL,
                     text.inside.symbol=NULL, text.beside.symbol=NULL,
                     text1.below.symbol=NULL, text2.below.symbol=NULL,
                     text3.below.symbol=NULL, text4.below.symbol=NULL,
                     filename=NULL, main=NULL, device="pdf", res=600, ...){
    device <- match.arg(device, c("pdf", "png", "plot"))
    ## device plot means we're just plotting, not generating a file.
    ## affected can be a vector or matrix with up to 4 columns.
    ## status (optional): 0 alive, missing, 1 dead.
    if(!is.null(is.adopted))
        warning("argument is.adopted not yet supported for kinship2 plotting.")
    if(!is.null(are.twins))
        warning("argument are.twins not yet supported for kinship2 plotting.")
    if(!is.null(are.consanguineous))
        warning("argument are.consanguineous not yet supported for kinship2 ",
                "plotting.")
    if(!is.null(is.sab.or.top))
        warning("argument is.sab.or.top not yet supported for kinship2 ",
                "plotting.")
    ## check the affected.
    if(!is.null(affected)){
        affected <- checkTransformNumeric(affected, name="affected")
        if(all(is.na(affected)))
            affected <- rep(0, length(individual))
    }else{
        affected <- rep(0, length(individual))
    }
    ## just creating the data.frame (although not using that one) to check all arguments.
    df <- buildHaplopaintDataframe(family=family, individual=individual,
                                   father=father, mother=mother, gender=gender,
                                   affected=affected, is.deceased=is.deceased,
                                   is.sab.or.top=is.sab.or.top,
                                   is.proband=is.proband, is.adopted=is.adopted,
                                   are.twins=are.twins,
                                   are.consanguineous=are.consanguineous,
                                   text.inside.symbol=text.inside.symbol,
                                   text.beside.symbol=text.beside.symbol,
                                   text1.below.symbol=text1.below.symbol,
                                   text2.below.symbol=text2.below.symbol,
                                   text3.below.symbol=text3.below.symbol,
                                   text4.below.symbol=text4.below.symbol
                                   )
    if(!is.null(is.deceased)){
        is.deceased <- checkTransformNumeric(is.deceased, name="is.deceased")
        ## should only have 0 and 1!
        is.deceased[is.na(is.deceased)] <- 0
        if(any(!(is.deceased) %in% c(0, 1)))
            stop("Only TRUE/FALSE or 1/0 allowed for is.deceased!")
    }else{
        is.deceased <- rep(0, length(individual))
    }
    ## OK, now trying to avoid a stupid bug in kinship2:
    if(any(is.na(affected)) & length(unique(affected)) == 2){
        ## append a dummy individual...
        affected <- c(affected, 0, 1)
        individual <- c(individual, "removeme1", "removeme2")
        father <- c(father, NA, NA)
        mother <- c(mother, NA, NA)
        gender <- c(gender, NA, NA)
        is.deceased <- c(is.deceased, 0, 0)
        ## have to add 2 elements to each other vector too.
        if(!is.null(text1.below.symbol))
            text1.below.symbol <- c(text1.below.symbol, NA, NA)
        if(!is.null(text2.below.symbol))
            text2.below.symbol <- c(text2.below.symbol, NA, NA)
        if(!is.null(text3.below.symbol))
            text3.below.symbol <- c(text3.below.symbol, NA, NA)
    }
    ## if we didn't get any error thus far, all should be right.
    ## generate a pedigree and plot, that's all.
    ped <- kinship2::pedigree(id=individual, dadid=father, momid=mother,
                              sex=sanitizeSex(gender))
    ## OK, let's start.
    FN <- filename
    if(device=="pdf"){
        if(is.null(FN))
            FN <- paste0(tempfile(), ".pdf")
        pdf(FN)
    }
    if(device=="png"){
        if(is.null(FN))
            FN <- paste0(tempfile(), ".png")
        png(FN)
    }
    ## checking additional input arguments.
    Args <- list(...)
    if(any(names(Args)=="cex")){
        Cex <- Args$cex
    }else{
        Cex <- 1
    }
    if(any(names(Args)=="symbolsize")){
        symbolsize <- Args$symbolsize
    }else{
        symbolsize <- 1
    }
    Coords <- plot(ped, status=is.deceased, affected=affected,
                   keep.par=TRUE, ...)
    oldxpd <- par()$xpd
    par(xpd=NA)
    ## if we do have any highlight ids or similar, i.e. text1.below.symbol etc.
    ## allow two: top left of symbol and top right of symbol.
    ## use the Coords$x, Coords$y, +/- Coords$boxw and pos=2, pos=4
    ## age below:
    ## use strheight to determine the size of the string.
    ## pos: y: Coords$y-Coores$boxh-strheight
    if(!is.null(text1.below.symbol)){
        strh <- strheight("0")
        text(x=Coords$x, y=(Coords$y+Coords$boxh-strh-strh/10), pos=1,
             labels=text1.below.symbol, cex=Cex-0.2)
    }
    if(!is.null(text2.below.symbol)){
        ## top left
        text(x=(Coords$x-Coords$boxw/8), y=Coords$y, pos=2,
             labels=text2.below.symbol, cex=Cex-0.2)
    }
    if(!is.null(text3.below.symbol)){
        ## top right
        text(x=(Coords$x+Coords$boxw/8), y=Coords$y, pos=4,
             labels=text3.below.symbol, cex=Cex-0.2)
    }
    if(!is.null(text4.below.symbol))
        warning("ks2paint does not support text4.below.symbol.")
    ## eventually, it might be nicer to highlight the id of probands in a
    ## different color... that's however not that easy as we see below: we do
    ## have to add the text for highlighted ids later to that plot to overplot
    ## the ids plotted in black
    if(!is.null(is.proband)){
        is.proband <- checkArgNumericLogical(is.proband, name="is.proband")
        Colors <- rep(NA, length(is.proband))
        Colors[is.proband==1] <- "#E41A1C"    ## RColorBrewer, Set1, (red)
        ## from kinship2: pretty tedious to get all the settings...
        psize <- par('pin')  # plot region in inches
        stemp1 <- strwidth("ABC", units='inches', cex=Cex)* 2.5/3
        stemp2 <- strheight('1g', units='inches', cex=Cex)
        stemp3 <- max(strheight(ped$id, units='inches', cex=Cex))
        xrange <- range(Coords$plist$pos[Coords$plist$nid >0])
        maxlev <- nrow(Coords$plist$pos)
        ht1 <- psize[2]/maxlev - (stemp3 + 1.5*stemp2)
        ht2 <- psize[2]/(maxlev + (maxlev-1)/2)
        wd2 <- .8*psize[1]/(.8 + diff(xrange))
        boxsize <- symbolsize* min(ht1, ht2, stemp1, wd2) # box size in inches
        vscale <- (psize[2]-(stemp3 + stemp2/2 + boxsize))/ max(1, maxlev-1)
        labh  <- stemp2/vscale   # height of a text string
        boxh <- Coords$boxh
        for(i in 1:maxlev){
            for(j in 1:Coords$plist$n[i]){
                k <- Coords$plist$nid[i,j]
                text(Coords$plist$pos[i,j], i + boxh + labh*.7, ped$id[k],
                     cex=Cex, adj=c(.5,1), col=Colors[k])
            }
        }
    }
    par(xpd=oldxpd)
    if(device!="plot")
        dev.off()
    ##return(Coords)
    invisible(FN)
}


#########################
## haplopaint uses the HaploPainter perl script/function.
haplopaint <- function(family=NULL, individual=NULL, father=NULL, mother=NULL,
                       gender=NULL, affected=NULL, is.deceased=NULL,
                       is.sab.or.top=NULL, is.proband=NULL, is.adopted=NULL,
                       are.twins=NULL, are.consanguineous=NULL,
                       text.inside.symbol=NULL, text.beside.symbol=NULL,
                       text1.below.symbol=NULL, text2.below.symbol=NULL,
                       text3.below.symbol=NULL, text4.below.symbol=NULL,
                       filename=NULL, main=NULL, device="pdf", res=600, ...){
    device <- match.arg(device, c("ps", "pdf", "svg", "png"))
    ## first check input arguments.
    df <- buildHaplopaintDataframe(family=family, individual=individual,
                                   father=father, mother=mother, gender=gender,
                                   affected=affected, is.deceased=is.deceased,
                                   is.sab.or.top=is.sab.or.top,
                                   is.proband=is.proband, is.adopted=is.adopted,
                                   are.twins=are.twins,
                                   are.consanguineous=are.consanguineous,
                                   text.inside.symbol=text.inside.symbol,
                                   text.beside.symbol=text.beside.symbol,
                                   text1.below.symbol=text1.below.symbol,
                                   text2.below.symbol=text2.below.symbol,
                                   text3.below.symbol=text3.below.symbol,
                                   text4.below.symbol=text4.below.symbol
                                   )
    rownames(df) <- df[, "individual"]
    ## handle childless founders:
    kins <- kinship(family=df[, "family"], id=df[, "individual"],
                    dadid=df[, "father"], momid=df[, "mother"])
    cols <- colSums(kins)
    if(any(cols <= 0.5)){
        kins <- kins[cols > 0.5, cols > 0.5]
        df <- df[colnames(kins), ]
        warning("Removed ", sum(cols <= 0.5), " childless founders!")
    }
    if(is.null(main))
        main <- df[1, "family"]
    if(is.null(filename))
        filename <- paste0(tempfile(), ".", device)
    ## save the data.frame to a temporary file
    dfFile <- tempfile()
    write.table(df, file=dfFile, sep="\t", row.names=FALSE, col.names=FALSE,
                quote=FALSE)
    ## call haplopaint
    plotcall <- paste0("perl ", options()$FamAgg$haplopaint, " -b -pedfile ",
                       dfFile, " -pedformat csv -outfile ", filename,
                       " -bgcolor \\#ffffff -outformat ", device,
                       " -resolution ", res, " -family ", main)
    res <- tryCatch(system(plotcall), error=function(e){return(e)})
    if(inherits(res, "simpleError")){
        stop("Error calling HaploPainter!", res$message)
    }else{
        if(res!=0)
            stop("Error calling HaploPainter! Please check error message.")
    }
    ## return the file name.
    invisible(filename)
}


#########################
## this function first checks all input parameters for correct input and
## subsequently builds a data.frame that can be used as input for haplopaint.
buildHaplopaintDataframe <- function(family=NULL, individual=NULL, father=NULL,
                                     mother=NULL, gender=NULL, affected=NULL,
                                     is.deceased=NULL, is.sab.or.top=NULL,
                                     is.proband=NULL, is.adopted=NULL,
                                     are.twins=NULL, are.consanguineous=NULL,
                                     text.inside.symbol=NULL,
                                     text.beside.symbol=NULL,
                                     text1.below.symbol=NULL,
                                     text2.below.symbol=NULL,
                                     text3.below.symbol=NULL,
                                     text4.below.symbol=NULL){
    if(is.null(individual) | is.null(father) |
       is.null(mother) | is.null(gender))
        stop("Arguments individual, father, mother and gender are required!")
    ninds <- length(individual)
    if(is.null(family))
        family <- rep(1, ninds)
    ## Fixing father & mother: Haplopaint requires 0
    father[is.na(father)] <- 0
    mother[is.na(mother)] <- 0
    ## now going through all arguments and, if present, adding them to the list
    argList <- list(family=family,
                    individual=as.character(individual),
                    father=as.character(father),
                    mother=as.character(mother)
                    )
    ## checking that sex has the correcto coding: 1:M, 2:F, 0:unknown
    gender <- as.numeric(sanitizeSex(gender))
    gender[is.na(gender)] <- 0
    argList <- c(argList, list(gender=gender))
    ## * affected: can be 0, 1 or any other number. each number is represented
    ##   by a different color: Note: we're adding +1 to each numeric affected
    ##   value! the coding is: NA -> 0, 0 -> 1, 1 -> 2!
    if(!is.null(affected)){
        affected <- checkTransformNumeric(affected, name="affected")
        if(length(which(affected < 0)))
            stop("Argument affected can not contain negative values!")
        if(any(is.na(affected))){
            ## adding 1 to each value and replacing NA with 0
            affected <- affected + 1
            affected[is.na(affected)] <- 0
        }else{
            ## adding +1 to each...
            affected <- affected + 1
        }
    }else{
        affected <- rep(0, ninds)
    }
    argList <- c(argList, list(affected=affected))
    ##
    ## * is.deceased
    ##   has to be either 0 or 1 (FALSE or TRUE)
    is.deceased <- checkArgNumericLogical(is.deceased, ninds,
                                          name="is.deceased")
    argList <- c(argList, list(is.deceased=is.deceased))
    ##
    ## * is.sab.or.top
    ##   has to be either 0 or 1 (FALSE or TRUE)
    is.sab.or.top <- checkArgNumericLogical(is.sab.or.top, ninds,
                                            name="is.sab.or.top")
    argList <- c(argList, list(is.sab.or.top=is.sab.or.top))
    ##
    ## * is.proband
    ##   has to be either 0 or 1 (FALSE or TRUE)
    is.proband <- checkArgNumericLogical(is.proband, ninds,
                                         name="is.proband")
    argList <- c(argList, list(is.proband=is.proband))
    ##
    ## * is.adopted
    ##   has to be either 0 or 1 (FALSE or TRUE)
    is.adopted <- checkArgNumericLogical(is.adopted, ninds,
                                         name="is.adopted")
    argList <- c(argList, list(is.adopted=is.adopted))
    ##
    ## * are.twins
    ##   NA, m_1, d_1 (m or d)_ followed by any text
    are.twins <- checkArgAreTwins(are.twins, ninds)
    argList <- c(argList, list(are.twins=are.twins))
    ##
    ## * are.consanguineous
    ##   NA and text specifying which individuals are consanguineous; always
    ##   two have to be consanguineous.
    are.consanguineous <- checkArgConsanguineous(are.consanguineous, ninds)
    argList <- c(argList, list(are.consanguineous=are.consanguineous))
    ##
    ## * text.inside.symbol
    ##   NA or any text
    if(!is.null(text.inside.symbol)){
        text.inside.symbol[is.na(text.inside.symbol)] <- ""
    }else{
        text.inside.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text.inside.symbol=text.inside.symbol))
    ##
    ## * text.beside.symbol
    ##   NA or any text
    if(!is.null(text.beside.symbol)){
        text.beside.symbol[is.na(text.beside.symbol)] <- ""
    }else{
        text.beside.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text.beside.symbol=text.beside.symbol))
    ##
    ## * text1.below.symbol
    ##   NA or any text
    if(!is.null(text1.below.symbol)){
        text1.below.symbol[is.na(text1.below.symbol)] <- ""
    }else{
        text1.below.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text1.below.symbol=text1.below.symbol))
    ##
    ## * text2.below.symbol
    ##   NA or any text
    if(!is.null(text2.below.symbol)){
        text2.below.symbol[is.na(text2.below.symbol)] <- ""
    }else{
        text2.below.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text2.below.symbol=text2.below.symbol))
    ##
    ## * text3.below.symbol
    ##   NA or any text
    if(!is.null(text3.below.symbol)){
        text3.below.symbol[is.na(text3.below.symbol)] <- ""
    }else{
        text3.below.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text3.below.symbol=text3.below.symbol))
    ##
    ## * text4.below.symbol
    ##   NA or any text
    if(!is.null(text4.below.symbol)){
        text4.below.symbol[is.na(text4.below.symbol)] <- ""
    }else{
        text4.below.symbol <- rep("", ninds)
    }
    argList <- c(argList, list(text4.below.symbol=text4.below.symbol))

    ## next we are going to evaluate whether we have the same length for all...
    ## if yes -> write temp data.frame and call haplopainter.
    Lengths <- unlist(lapply(argList, length))
    if(length(unique(Lengths)) != 1){
        cat("Got lengths: ", paste(paste(names(argList), ": ", Lengths),
                                   collapse=", "),  "\n")
        stop("All arguments have to have the same number of elements!")
    }
    do.call(cbind, argList)
}

## tries to transform x to numeric and, in case it's not possible, throws
## an error.
checkTransformNumeric <- function(x, name=""){
    tryCatch(x <- as.numeric(x),
             error=function(e){
                 stop("Error while checking argument ", name,
                      ": argument is not numeric!")
             },
             warning=function(w){
                 stop("Error while checking argument ", name,
                      ": argument is not numeric!")
             })
    x
}

## check the are.twins argument.
checkArgAreTwins <- function(are.twins=NULL, ninds){
    if(!is.null(are.twins)){
        if(is.character(are.twins))
            are.twins <- factor(are.twins)
        if(!is.factor(are.twins))
            stop("are.twins has to be either a character vector or factor!")
        ## checking the levels: have to start with m_ or d_
        twinlevels <- levels(are.twins)
        if(!all(substring(twinlevels, 1, 2) %in% c("m_", "d_")))
            stop("elements in are.twins have to be either NA or have to ",
                 "start with m_ (monozygotic twins) or d_ (dizygotic twins)!")
        ## each level has to be present at least twice.
        if(any(table(are.twins) < 2))
            stop("each text to specify twins has to be present at least ",
                 "twice in argument are.twins!")
        ## OK, seems to be OK
        are.twins <- as.character(are.twins)
        are.twins[is.na(are.twins)] <- ""
    }else{
        are.twins <- rep("", ninds)
    }
    are.twins
}

##
## checks whether the input argument can be transformed into a numeric with
## values 0 and 1.
checkArgNumericLogical <- function(x=NULL, ninds, name=""){
    if(!is.null(x)){
        x <- checkTransformNumeric(x, name=name)
        if(any(!(x %in% c(0, 1))))
            stop(name, " should be a logical vector or a numeric vector ",
                 "with 0 (no) and 1 (yes)!")
    }else{
        x <- rep(0, ninds)
    }
    x
}

##
##
##
checkArgConsanguineous <- function(x=NULL, ninds){
    if(!is.null(x)){
        x <- factor(x)
        ## each level has to be present exactly twice.
        if(any(table(x) != 2)){
            stop("Each text specifying consanguineous couples in argument ",
                 "are.consanguineous should be present exactly twice!")
        }
        x <- as.character(x)
        x[is.na(x)] <- ""
    }else{
        x <- rep("", ninds)
    }
    x
}
