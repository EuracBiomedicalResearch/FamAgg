####============================================================
##  import ped files
##
##  Function to import the ped file format used in plink: Expected
##  columns:
##  + Family ID
##  + Sample ID
##  + Paternal ID (0 if father isn't in dataset)
##  + Maternal ID (0 if mother isn't in dataset)
##  + Sex (1=male; 2=female; other=unknown)
##  + Affection (0=unknown; 1=unaffected; 2=affected)
##  + Genotypes (space or tab separated, 2 for each marker. 0, -9,
##    non-numeric=missing)
##  We will ignore the Genotype information though.
####------------------------------------------------------------
doImportPed <- function(x, sep="\t", header=FALSE, ...){
    df <- read.table(x, sep=sep, header=header)
    doProcessPed2Df(df)
}

####============================================================
##  import fam files
##
##  Function to import the fam file format used in plink. The fam
##  file corresponds to the first 6 columns of the ped file.
####------------------------------------------------------------
doImportFam <- function(x, sep="\t", header=FALSE, ...){
    df <- read.table(x, sep=sep, header=header)
    doProcessPed2Df(df)
}

####============================================================
##  import generic txt files
##
####------------------------------------------------------------
doImportGeneric <- function(x, sep="\t", header=TRUE, id.col="id",
                            family.col="family", father.col="father",
                            mother.col="mother", sex.col="sex", ...) {
    df <- read.table(x, sep=sep, header=header)
    if (header) {
        ## check if I have all required column names...
        if(!all(c( id.col, family.col, father.col, mother.col, sex.col) %in%
                colnames(df)))
            stop("One or more required column names not found! Expect: ",
                 family.col, ", ", id.col, ", ", father.col, ", ",
                 mother.col, ", ", sex.col, "!")
        df <- df[, c(family.col, id.col, father.col, mother.col, sex.col)]
        colnames(df) <- .PEDCN
    }
    df
}

####============================================================
##  process a data.frame corresponding to a ped file and re-format
##  it according to our needs. The expected columns are:
##  + Family ID
##  + Sample ID
##  + Paternal ID (0 if father isn't in dataset)
##  + Maternal ID (0 if mother isn't in dataset)
##  + Sex (1=male; 2=female; other=unknown)
##  + Affection (0, -9, non-numeric=unknown; 1=unaffected; 2=affected)
##
####------------------------------------------------------------
doProcessPed2Df <- function(x){
    ## Subset the ped data.frame to it's first 6 columns.
    if(ncol(x) < 6)
        stop("A fam/ped file is expected to have at least 6 columns! ",
             "The present one has however only ", ncol(x), "!")
    ## Extract family ID.
    FamId <- as.character(x[, 1])
    ## Extract individual ID.
    IId <- as.character(x[, 2])
    ## Checking father id and mother id; should be ALL in the pedigree!
    FatherId <- as.character(x[, 3])
    ## Set all empty ("") and 0s to NA
    FatherId[which(FatherId == "" | FatherId == "0" )] <- NA
    if(!all(FatherId[!is.na(FatherId)] %in% IId))
        stop("If provided (i.e. different from '0') the father ID has",
             " to correspond to an ID in the pedigree!")
    MotherId <- as.character(x[, 4])
    MotherId[which(MotherId == "" | MotherId == "0" )] <- NA
    if(any(!(MotherId[!is.na(MotherId)] %in% IId)))
        stop("If provided (i.e. different from '0') the mother ID has",
             " to correspond to an ID in the pedigree!")
    ## Check the sex: 1=male, 2=female; other=unknown.
    Sex <- rep(NA, nrow(x))
    Sex[which(x[, 5] == "1")] <- 1
    Sex[which(x[, 5] == "2")] <- 2
    ## Affected:
    ## 0, -9 or any non-numeric are mapped to NA (unknown)
    ## 1 is mapped to 0 (unaffected)
    ## 2 is mapped to 1 (affected)
    Affected <- x[, 6]
    ## Check if we've got case/control or numeric. We're not supporting the latter!
    if(is.numeric(Affected)){
        amNa <- is.na(Affected)
        if(!all(Affected[!amNa] %in% c(1, 2, 0, -9)))
            stop("FamAgg supports only categorical case/control phenotype",
                 " information (i.e. values 1, 2, 0 or -9 in phenotype ",
                 "column 6 of the data file)!")
    }
    ## Convert to character and fix to only NA, 0 and 1.
    Affected <- as.character(Affected)
    NewAff <- rep(NA, length(Affected))
    NewAff[Affected == "1"] <- 0
    NewAff[Affected == "2"] <- 1
    ## Compile the data.frame
    df <- data.frame(family=FamId, id=IId, father=FatherId, mother=MotherId,
                     sex=Sex, trait=NewAff, stringsAsFactors=FALSE)
    df$sex <- sanitizeSex(df$sex)
    df
}

####============================================================
##  Format a pedigree data.frame as extracted from FAData in "ped"
##  format.
####------------------------------------------------------------
doProcessDf2Ped <- function(x){
    ## Assuming correct column names (family, id, father, mother, sex,
    ## affected)
    if(!any(colnames(x) == "affected"))
        x <- cbind(x, affected=rep(NA, nrow(x)))
    x <- x[, c("family", "id", "father", "mother", "sex", "affected")]
    x$father[is.na(x$father)] <- 0
    x$mother[is.na(x$mother)] <- 0
    x$sex <- as.numeric(x$sex)
    x$sex[is.na(x$sex)] <- 0
    x$affected <- x$affected + 1
    x$affected[is.na(x$affected)] <- 0
    x
}


####============================================================
##  doImport
##
##  The main import function that switches based on the file extension.
####------------------------------------------------------------
doImport <- function(file, ...){
    extension <- .fileExtension(file)
    file <- path.expand(file)
    if(extension == "fam")
        return(doImportFam(file, ...))
    if(extension == "ped")
        return(doImportPed(file, ...))
    ## Fall back to "generic import"
    doImportGeneric(file, ...)
}

.fileExtension <- function(file){
    if(!grepl("\\.", file))
        stop("Unable to identify extension for file '", file, "'")
     ext <- sub(".*\\.", "", sub("\\.gz$|\\.gzip$", "", basename(file)))
    if(ext=="")
        stop("Unable to identify extension for file '", file, "'")
    tolower(ext)
}



