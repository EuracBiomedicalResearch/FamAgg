\name{FAData-class}
\Rdversion{1.1}
\docType{class}
\alias{FAData-class}
\alias{FAData}
\alias{$}
\alias{$,FAData-method}
\alias{[,FAData,ANY-method}
\alias{[,FAData,ANY,ANY,ANY-method}
\alias{affectedIndividuals}
\alias{affectedIndividuals,FAData-method}
\alias{age}
\alias{age,FAData-method}
\alias{age<-}
\alias{age<-,FAData-method}
\alias{buildPed}
\alias{buildPed,FAData-method}
\alias{export}
\alias{export,FAData-method}
\alias{family}
\alias{family,FAData-method}
\alias{kinship}
\alias{kinship,FAData-method}
\alias{pedigree}
\alias{pedigree,FAData-method}
\alias{pedigree,missing-method}
\alias{pedigree<-}
\alias{pedigree<-,FAData-method}
\alias{pedigreeSize}
\alias{pedigreeSize,FAData-method}
\alias{phenotypedIndividuals}
\alias{phenotypedIndividuals,FAData-method}
\alias{plotPed}
\alias{plotPed,FAData-method}
\alias{show}
\alias{show,FAData-method}
\alias{trait}
\alias{trait,FAData-method}
\alias{trait<-}
\alias{trait<-,FAData-method}

\title{Pedigree data information}

\description{
  \code{FAData} objects conveniently store pedigree along with trait
  information. This object is the central data structure from the
  \code{FamAgg} package. Basic usage pedigree analysis methods are
  described on this page and on the \code{\link{PedigreeUtils}} help
  page, familial aggregation analysis methods on the
  \code{\link{PedigreeAnalysis}} help page.

  See the section about the \code{pedigree} \code{data.frame} below for
  a detailed description of the encoding of missing trait data or
  founder individuals in \code{FamAgg}.
}

\section{Objects from the Class}{
  \code{FAData} objects are created by the constructor function
  \code{FAData} and should not be directly created by a call to
  \code{new}.
}
\usage{

\S4method{affectedIndividuals}{FAData}(object)

\S4method{age}{FAData}(object)

\S4method{age}{FAData}(object) <- value

\S4method{buildPed}{FAData}(object, id=NULL, family = NULL, max.generations.up=3,
                            max.generations.down=16, prune=FALSE, ...)

\S4method{export}{FAData}(object, con, format="ped", ...)

FAData(pedigree, age, trait, traitName, header=FALSE, sep="\t", id.col="id",
       family.col="family", father.col="father", mother.col="mother",
       sex.col="sex")

\S4method{family}{FAData}(object, id=NULL, family=NULL,
                          return.type="data.frame")

\S4method{kinship}{FAData}(id, ...)

\S4method{pedigree}{FAData}(object, return.type="data.frame")

\S4method{pedigree}{FAData}(object) <- value

\S4method{pedigreeSize}{FAData}(object)

\S4method{phenotypedIndividuals}{FAData}(object)

\S4method{plotPed}{FAData}(object, id=NULL, family=NULL, filename=NULL,
                           device="plot", symbol.related=NA,
                           proband.id=NULL, highlight.ids=NULL,
                           only.phenotyped=FALSE,
                           label1=age(object), label2=NULL, label3=NULL,
                           ...)

\S4method{show}{FAData}(object)

\S4method{trait}{FAData}(object, na.rm=FALSE)

\S4method{trait}{FAData}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{age}{
    For \code{FAData}: either a character(1) specifying the file name
    from which the age should be read or a named numeric vector of ages
    with the names corresponding to the ids of the individuals in the
    pedigree.
  }

  \item{con}{
    For \code{export}: the file name or connection to a file to which
    the pedigree information should be exported.
  }

  \item{device}{
    For \code{plotPed}: the device of file format in which the plot
    should be saved. See details for allowed values.
  }

  \item{family}{
    For \code{buildPed}: the id of the family for which the pedigree
    should be returned.
    For \code{family}: the id of the family for which the pedigree
    should be returned (full pedigree of the family).
    For \code{plotPed}: the id of the family for which the pedigree
    should be plotted.
  }

  \item{family.col}{
    For \code{FAData}: the name of the column containing the id of the
    families, if argument \code{pedigree} is a \code{character(1)}
    specifying the file name from which the pedigree should be read and
    \code{header} is \code{FALSE}.
  }

  \item{father.col}{
    For \code{FAData}: the name of the column containing the id of the
    father, if argument \code{pedigree} is a \code{character(1)}
    specifying the file name from which the pedigree should be read and
    \code{header} is \code{FALSE}.
  }

  \item{filename}{
    For \code{plotPed}: a character string specifying the name of the
    file to which the plot should be saved. If none is submitted, the
    plot is saved to a temporary file.
  }

  \item{format}{
    For \code{export}: the format in which the pedigree should be
    exported. At present only \code{"ped"} and \code{"fam"} are
    exported, i.e. the file formats from plink
    (\url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml}).
  }

  \item{header}{
    For \code{FAData}: only used if argument \code{pedigree} is a
    \code{character(1)}, i.e. the file name from which the pedigree
    should be read. The \code{header} argument is passed to the
    \code{read.table} function, i.e. should be set to \code{TRUE} if the
    file contains column headers.
  }

  \item{highlight.ids}{
    A list of character vector(s) of ids that should be labeled. The
    name(s) of the character vector(s) is/are used as the text to label
    the individuals (the text is shown below the symbol of the
    individuals). Up to 3 character vectors are
    supported. Alternatively, a single character vector of ids can be
    submitted in which case the individuals are labeled with an asterisc
    (\code{"*"}).
  }

  \item{id}{
    For method \code{kinship}: the \code{FAData} object from which
    the kinship matrix should be extracted, for all other methods the id
    of the individual.

    For method \code{plotPed}: the id of the individual for which the
    pedigree should be built (see \code{buildPed}) and plotted.

    Note: \code{id} can be a \code{numeric} or a
    \code{character}. Numeric ids are internally converted to
    \code{character}.
  }

  \item{id.col}{
    For \code{FAData}: the name of the column containing the id of the
    individuals, if argument \code{pedigree} is a \code{character(1)}
    specifying the file name from which the pedigree should be read and
    \code{header} is \code{FALSE}.
  }

  \item{label1}{
    For \code{plotPed}: labels that should be plotted below the symbol
    for each individual. Should be either a named vector with names
    corresponding to the ids of the individuals in the pedigree or a
    vector of the same length than individuals that are to be
    plotted. For the former it is sufficient to just specify the labels
    for the individuals that should be shown.
  }

  \item{label2}{
    For \code{plotPed}: see \code{label1}. The labels are plotted in the
    second line below the symbol if \code{HaploPainter} is used to
    generate the plot, or on the top left corner of the individual's
    symbol for \code{kinship2} plotting.
  }

  \item{label3}{
    For \code{plotPed}: see \code{label1}. The labels are plotted in the
    third line below the symbol if \code{HaploPainter} is used to
    generate the plot, or on the top right corner of the individual's
    symbol for \code{kinship2} plotting.
  }

  \item{max.generations.down}{
    For \code{buildPed}: the maximal number of generations to look for
    children.
  }

  \item{max.generations.up}{
    For \code{buildPed}: the maximal number of generations to look for
    ancestors.
  }

  \item{mother.col}{
    For \code{FAData}: the name of the column containing the id of the
    mother, if argument \code{pedigree} is a \code{character(1)}
    specifying the file name from which the pedigree should be read and
    \code{header} is \code{FALSE}.
  }

  \item{na.rm}{
    For \code{trait}: whether missing values in trait should be returned
    or not.
  }

  \item{object}{
    The \code{FAData} object.
  }

  \item{only.phenotyped}{
    Wheter only phenotyped individuals, i.e. individuals with a non-NA
    value in column \code{affected} (the trait information). Requires
    this information to be present.
  }

  \item{pedigree}{
    For \code{FAData}: either a \code{data.frame} with the pedigree
    information or a \code{character(1)} specifying the file name from
    which the pedigree should be read. See description below for more
    details.
  }

  \item{proband.id}{
    For \code{plotPed}: character vector with the id(s) of one ore more
    individuals that should be highlighted as probands.
    \code{HaploPainter} indicates probands with a \code{"P"} next to the
    symbol and an arrow pointing to the symbol.
  }

  \item{prune}{
    For \code{buildPed}: whether the smallest possible (connected)
    pedigree for the submitted ids should be build. This makes only
    sense if more than one id is submitted.
  }

  \item{return.type}{
    Either \code{"data.frame"} or \code{"pedigree"} if the pedigree
    information should be returned as a \code{data.frame} or
    \code{pedigreeList} object as defined in the \code{kinship2}
    package.
  }

  \item{sep}{
    For \code{FAData}: only used if argument \code{pedigree} is a
    \code{character(1)}, i.e. the file name from which the pedigree
    should be read. The \code{sep} argument is passed to the
    \code{read.table} function and specifies the field separator.
  }

  \item{sex.col}{
    For \code{FAData}: the name of the column spefifying the sex of the
    individuals, if argument \code{pedigree} is a \code{character(1)}
    representing the file name from which the pedigree should be read and
    \code{header} is \code{FALSE}.
  }

  \item{symbol.related}{
    For \code{plotPed}: the symbol which should be used to label
    individuals sharing kinship with the id for which the pedigree is
    generated and plotted.
  }

  \item{trait}{
    For \code{FAData}: a numeric vector with \code{0}, \code{1} and
    \code{NA} or a logical vector indicating unaffected (but
    phenotyped), affected and not phenotyped individuals.
  }

  \item{traitName}{
    For \code{FAData}: an optional name for the trait.
  }

  \item{value}{
    For \code{age<-}: a named numeric vector. The names (at least some
    of them) have to match the ids in the pedigree of the object.

    For \code{pedigree<-}:

    For \code{trait<-}: a named numeric vector with \code{0}, \code{1}
    and \code{NA} or a logical vector with \code{FALSE}, \code{TRUE},
    \code{NA} for not affected, affected and not tested. The names (at
    least some of them) have to match the ids in the pedigree of the
    object.

  }

  \item{...}{
    Additional arguments to be passed to the functions
    (\code{\link{doPlotPed}} for \code{plotPed}).
  }

}
\section{Slots}{
  \describe{
    \item{age}{
      A (named) numerical vector with the age of the individuals. It is
      suggested to use the getter and setter methods described below to
      access this slot.
    }

    \item{pedigree}{
      A \code{data.frame} with the pedigree. It is suggested to use the
      getter and setter methods described below to access this slot.
    }

    \item{.kinship}{
      The kinship matrix for the kinship of each individual in the
      pedigree with each other. This slot should not be accessed
      directly, but the \code{kinship} method should be used
      instead.
    }

    \item{traitname}{
      The name of the trait being stored in the object.
    }

    \item{.trait}{
      A numerical vector with the trait information, \code{0}, \code{1},
      \code{NA}, for phenotyped but not affected, affected and not
      tested, respectively. This slot should not be accessed directly,
      but the \code{trait} and \code{trait<-} methods should be used
      instead that ensure that the data is matched to the information in
      the pedigree.
    }

  }
}

\section{Constructors, importing and exporting data}{
  \describe{
    \item{FAData}{
      Constructor function to create a new \code{FAData} instance. In
      addition to submitting the pedigree information as
      \code{data.frame}, \code{pedigree} or \code{pedigreeLinst} it is
      possible to specify the name of the file from which the pedigree
      information should be read.
      The recognizes and imports plink \emph{ped} and \code{fam} files
      (\url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml}) or from
      generic text files.
      For the latter, arguments
      \code{header}, \code{sep}, \code{family.col}, \code{id.col},
      \code{father.col}, \code{mother.col} and \code{sex.col} allow to
      further specify which columns of the file contain what information
      etc.
      If argument \code{pedigree} is a \code{data.frame}, the column
      names \code{"family"}, \code{"id"}, \code{"father"},
      \code{"mother"} and \code{"sex"} are expected. Any additional
      columns are dropped.

      The sex is expected to be encoded either as a numeric 1 (male),
      2 (female) with any other number or \code{NA} representing
      unknown, or as a character vector or factor with \code{"M"},
      \code{"m"}, \code{"Male"} or \code{"male"} for male and \code{"F"},
      \code{"f"}, \code{"Female"} or \code{"female"} for female.
    }

    \item{export}{
      Export pedigree data to a file.
    }
  }
}

\section{Accessors and subsetting}{
  \describe{

    \item{object$name}{
      Access \code{name} column in the pedigree of the \code{FAData}
      object. The function returns a named vector wirh the names
      corresponding to the ids of the individuals or \code{NULL} if
      \code{name} does not correspoind to a column name in the pedigree.
      The trait data can be accessed either by \code{object$trait} or
      \code{object$affected}.
    }

    \item{age}{
      Returns the age of the individuals as a named numeric vector.
      If the pedigree is set, the order of the values corresponds always
      to the ordering of the individuals in the pedigree with \code{NA}
      for individuals for which the age is unknown. In case the age was
      never set it returns a vector of NAs with length equal to the
      number of individuals.
    }

    \item{age<-}{
      Setter for the age. Value has to be a named numeric vector.
    }

    \item{pedigree}{
      Returns the pedigree either as a \code{data.frame} or a
      \code{pedigreeList} object (defined in the \code{kinship2}
      package) depending on the value of the parameter
      \code{return.type} (i.e. either \code{return.type="data.frame"} or
      \code{return.type="pedigree"}). If \code{pedigree} is called on
      any other object than a \code{FAData} object (or any object
      that inherits from that object), the
      \code{\link[kinship2]{pedigree}} method from the
      \code{kinship2} package is called.

      For the default return type (i.e. \code{return.type="data.frame"})
      a \code{data.frame} is returned with the following columns:
      \code{"family"}: the ID of the family,
      \code{"id"}: the ID of the individual,
      \code{"father"}: the ID of the individual's father. Founder
      individuals, i.e. individuals for whom the father and mother is
      not known in the data set, contain a \code{NA} in this column.
      \code{"mother"}: the ID of the individual's father. Founder
      individuals, i.e. individuals for whom the father and mother is
      not known in the data set, contain a \code{NA} in this column.
      \code{"sex"}: the sex of the individuals encoded as a
      \code{factor} with levels \code{"M"} and \code{"F"} for male and
      female, or \code{NA} for not known.
      If trait information is available in the object the returned
      \code{data.frame} will also contain a column named \code{affected}
      with the information whether the individual is affected
      (\code{1}), not affected (\code{0}) or was not tested/phenotyped
      \code{NA}.
    }

    \item{pedigree<-}{
      Setter for the pedigree slow. Value can be a \code{data.frame}
      with columns containing the family id, individual id, father id,
      mother id and sex (in this order) or a \code{pedigree} or
      \code{pedigreeList} object as defined in the \code{kinship2}
      package.
    }

    \item{object[i, ]}{
      Subsets the \code{FAData} object to individuals specified with
      \code{i} which can be a logical, numeric or character vector. For
      the latter, the elements have to be the ids of the individuals
      (i.e. rownames of \code{pedigree(object)}). Returns the sub-setted
      object.
      Note that subsetting other than by family might result in a
      non-valid pedigree (e.g. if mother or father ID are not available
      in the sub-setted pedigree).
    }

    \item{trait}{
      Get the trait vector from the object. By default, the ordering
      is the same as pedigree, setting argument \code{na.rm=TRUE}
      removes all \code{NA} values, thus the ordering and length might
      be different. Returns a named vector with the names corresponding
      to the ids of the individuals.
    }

    \item{trait<-}{
      Setter for the trait slot. Can be a named numeric vector (values
      \code{0}, \code{1} and \code{NA}) or logical vector (values
      \code{FALSE}, \code{TRUE} and \code{NA}) with the names matching
      the ids of the individuals in the pedigree. The method internally
      matches and re-orders the trait vector to match the ordering of
      the ids in the pedigree.
    }

  }
}

\section{Basic usage}{
  \describe{

    \item{affectedIndividuals}{
      Returns a character vector with the ids of the affected
      individuals, i.e. the id of the individuals with a value other
      than 0 or \code{NA} in the trait.
      If no trait data is available the method returns \code{NULL}.
    }

    \item{buildPed}{
      Builds a pedigree for the specified id(s) containing generations
      defined by \code{max.generations.up} and
      \code{max.generations.down} and returns it as a
      \code{data.frame}. The pedigree contains all individuals in the
      family sharing kinship with the input individual(s) and mates
      needed to complete the pedigree.
      For \code{prune=TRUE} the function tries to
      find the smallest connected pedigree for all the submitted ids.
    }

    \item{family}{
      Returns the pedigree for a full family. In contrast to
      \code{buildPed} which constructs a (sub)pedigree for a specific
      individual, this method returns the pedigree of the complete
      family for an individual (if \code{id} is specified). The function
      returns either a \code{data.frame} or a \code{pedigreeList} with
      the pedigree for the family.
    }

    \item{kinship}{
      Extracts the pre-calculated kinship matrix, i.e. a symmetric
      matrix with the kinship between all individuals in the
      pedigree. The matrix is calculated using the \code{kinship} method
      provided by the \code{kinship2} package [Sinwell (2014)].
      The function returns a \code{dsCMatrix} from the \code{Matrix}
      package.
    }

    \item{pedigreeSize}{
      Returns the size, i.e. the number of individuals (rows) in the
      pedigree.
    }

    \item{phenotypedIndividuals}{
      Returns a character vector with the ids of the phenotyped
      individuals, i.e. the id of all individuals that have a non-NA
      value in thetrait.
      If no trait data is available the method returns \code{NULL}.
    }

    \item{plotPed}{
      Creates the pedigree for the submitted id(s) or family and plots
      it (i.e. saves it to the specified file). See details above for
      more information.
      Returns the file name of the file to which the pedigree plot was
      exported or \code{NULL} for \code{kinship2} plotting and
      \code{device="plot"}. See \code{\link{doPlotPed}} for more
      information.
    }

  }
}
\section{Pedigree analysis methods}{
  Methods for familial aggregation and other pedigree analysis methods
  are described on the \link{PedigreeAnalysis} help page.
}

\section{Pedigree utilities}{
  A variety of different pedigree utilities are defined for
  \code{FAData} objects. For the full list of methods see the
  \link{PedigreeUtils} help page.
}

\details{

  See sections below for a description of the individual methods.

  The \code{buildPed} method is a combination of the methods
  \code{getAncestors}, \code{getChildren} and \code{getMissingMate},
  i.e. it first gets all ancestors for the specified id(s), determines
  then the children of all of the ids (submitted ids and their
  ancestors) and at last looks for any missing mates/spouses to complete
  the pedigree.

  The \code{plotPed} function uses either the external perl program
  \code{HaploPainter} or the plotting capabilities of the
  \code{kinship2} package.
  With \code{HaploPainter}, as it is an external too, it is not possible
  to display the plot directly, but each plot is automatically saved to
  a file (either \code{"pdf"}, \code{"ps"}, \code{"svg"} or
  \code{"png"}; can be specified with the \code{device} parameter).
  Plotting with \code{kinship2} (the default) allows to display the plot
  (\code{device="plot"}) or export it to a file (\code{device="pdf"} or
  \code{device="png"}).

  The \code{\link{switchPlotfun}} function can be used to change the plotting
  system.

}
\note{
  The ids of individuals, father, mother and family can be either
  numeric or characters, internally, all ids will however be handled as
  characters.

  The \code{pedigree<-} setter method removes all white spaces in
  columns \code{"id"}, \code{"family"}, \code{"father"} and
  \code{"mother"} of the pedigree.
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer.
}
\references{
  Sinwell JP, Therneau TM & Schaid DJ (2014) The kinship2 R package for
  pedigree data. \emph{Human heredity} 78:91-93.
}
\seealso{
  \code{\link[kinship2]{pedigree}},
  \code{\linkS4class{FAProbResults}},
  \code{\linkS4class{FAKinGroupResults}},
  \code{\linkS4class{FAKinSumResults}},
  \code{\linkS4class{FAGenIndexResults}},
  \code{\link{doPlotPed}},
  \code{\link{PedigreeUtils}},
  \code{\link{getAll}},
  \code{\link{PedigreeAnalysis}}
}
\examples{
##########################
##
##  Create a new FAData object
##
## Load the Minnesota Breast Cancer record and subset to the
## first families.
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid==4 | minnbreast$famid==5, ]
mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
## Renaming column names
colnames(mbped) <- c("family", "id", "father", "mother", "sex")
## Defining the optional argument age.
Age <- mbsub$endage
names(Age) <- mbsub$id
## Create the object
fad <- FAData(pedigree=mbped, age=Age)

fad

## Extract the ids directly...
head(fad$id)

## Extract the kinship matrix
dim(kinship(fad))

## What's the size of the pedigree?
pedigreeSize(fad)

## Importing a "ped" file.
pedFile <- system.file("txt/minnbreastsub.ped.gz", package="FamAgg")
## Quick glance at the file.
readLines(pedFile, n=1)
fad <- FAData(pedFile)

head(pedigree(fad))

## Creating the FAData reading data from a txt file.
pedFile <- system.file("txt/minnbreastsub.txt", package="FamAgg")
fad <- FAData(pedigree=pedFile, header=TRUE, id.col="id",
              family.col="famid", father.col="fatherid",
              mother.col="motherid")
## Adding the age
age(fad) <- Age
fad
## List all families in the pedigree along with the number of
## individuals
table(fad$family)

##########################
##
##  Basic usage
##
## Extracting the pedigree information
ped <- pedigree(fad)
## By default the pedigree is returned as a data.frame.
class(ped)
head(ped)

## In addition, we can extract the pedigree as a pedigreeList
pedigree(fad, return.type="pedigree")

## Return the ids of all ancestors of individual 6
## up to 3 generations
getAncestors(fad, id="6")

## Build the pedigree for individual 6: this includes all of its
## children and all of its ancestors up to the maximal number of
## specified generations.
buildPed(fad, id=6)
## Which is a sub-pedigree of the complete family:
family(fad, id=6)

## In addition we can specify manually some ids in the pedigree and
## generate the smallest possible pedigree containing all ids:
buildPed(fad, id=c(6, 23, 28), prune=TRUE)

## Get the list of all ids sharing kinship with individuals
## 5 and 9
shareKinship(fad, id=c("5", "9"))

## Subset the fad to family "4"
subFad <- fad[fad$family == "4", ]
subFad

## Export the pedigree from this family to a ped file
tmpFile <- tempfile()
export(subFad, con=tmpFile, format="ped")

head(read.table(tmpFile, sep="\t"))

##########################
##
##  Plotting
##
## Plot the pedigree for individual 6.
plotPed(fad, id=6)

## Alternatively, exporte it to a temporary file
pfile <- plotPed(fad, id=6, device="pdf")
pfile

## Highlighting some of the individuals:
## first get to know which other individuals are in the pedigree
plotPed(fad, id=6, highlight.ids=list(hello=c(1, 2, 4)))


##########################
##
##  Adding trait data
##
fad <- FAData(pedigree=mbped, age=Age)
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id
trait(fad) <- tcancer
## Now we can plot the pedigree also showing the affected status.
plotPed(fad, id=6)


## Alternatively, create the FAData with the trait data
fad <- FAData(pedigree=mbped, trait=mbsub$cancer, traitName="cancer")
plotPed(fad, id=6)

}
\keyword{classes}

