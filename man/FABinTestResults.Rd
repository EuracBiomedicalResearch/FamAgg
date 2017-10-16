\name{FABinTestResults-class}
\Rdversion{1.1}
\docType{class}
\alias{FABinTestResults-class}
\alias{FABinTestResults}
\alias{[,FABinTestResults,ANY-method}
\alias{[,FABinTestResults,ANY,ANY,ANY-method}
\alias{result,FABinTestResults-method}
\alias{trait<-,FABinTestResults-method}


\title{Binomial test for familial aggregation}

\description{
  The \code{FABinTestResults} object contains the results from a
  simple binomial to test whether the number of affected in a trait are
  higher than expected by chance. For more details on the method please
  see \code{\link{binomialTest}}.

}
\section{Objects from the Class}{
  \code{FABinTestResults} objects are returned by the
  \code{\link{binomialTest}} function.
}
\usage{

\S4method{result}{FABinTestResults}(object, method="BH")

\S4method{trait}{FABinTestResults}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{method}{
    The multiple hypothesis testing method. All methods supported by
    \code{\link{p.adjust}} are allowed.
  }

  \item{object}{
    The \code{FABinTestResults} object.
  }

  \item{value}{
    For \code{trait<-}: can be a named numeric, character or factor
    vector. The names (at least some of them) have to match the ids in
    the pedigree of the object.
  }

}
\section{Extends}{
  Class \code{\linkS4class{FAData}} directly.
}
\section{Slots}{
  \describe{
    \item{result}{
      The results \code{data.frame}.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{result}{
      Returns the result from the test as a \code{data.frame} with
      columns:
      
      \code{"trait_name"}: the name of the trait.
      
      \code{"total_phenotyped"}: total number of phenotyped individuals in
      the trait.
      
      \code{"total_affected"}: total number of affected individuals in the
      trait.
      
      \code{"family"}: the family id. If a global test is used (i.e. if
      the pedigree consists of a single family, or \code{global = TRUE}
      was provided, the column shows \code{"full pedigree"}).
      
      \code{phenotyped}: the number of phenotyped individuals in the
      family.

      \code{affected}: the number of affected individuals in the family.

      \code{pvalue}: the p-value from the binomial test (conducted using
      the \code{\link{binom.test}} function).

      \code{prob}: the probability of being affected. Either a
      \emph{local} probability calculated based on all affected and
      phenotyped individuals in the whole pedigree, or a \emph{global}
      (population) probability that has to be provided with argument
      \code{prob}.

      \code{padj}: the p-value adjusted for multiple hypothesis testing
      using the method defined with argument \code{"method"}.
    
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\details{

  A call to the setter methods \code{trait<-} resets any test
  results present in the \code{result} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.

}
\note{
  Subsetting (using the \code{[} operator) is not supported.
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer, Christian Weichenberger
}
\seealso{
  \code{\linkS4class{FAData}},
  \code{\link{kinship}},
  \code{\link{trait}},
  \code{\link{probabilityTest}},
  \code{\link{kinshipGroupTest}},
  \code{\link{kinshipSumTest}},
  \code{\link{genealogicalIndexTest}},
  \code{\link{familialIncidenceRateTest}},
  \code{\link{fsirTest}},
  \code{\link{plotPed}}
}
\examples{
##########################
##
##  Perform the analysis
##
## Load the test data.
data(minnbreast)

## Subset to some families and generate a pedigree data.frame.
mbsub <- minnbreast[minnbreast$famid == 4 | minnbreast$famid == 5 |
                    minnbreast$famid == 6 | minnbreast$famid == 7 |
                    minnbreast$famid == 8, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- c("family", "id", "father", "mother", "sex")

## Generate the FAData.
fad <- FAData(pedigree=PedDf)

## Specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## Perform the test:
bir <- binomialTest(fad, trait=tcancer, traitName="cancer")
result(bir)

## Calculating the probability of being affected from the whole data set.
prob <- sum(minnbreast$cancer, na.rm = TRUE) / sum(!is.na(minnbreast$cancer))
bir <- binomialTest(fad, trait = tcancer, prob = prob)
result(bir)

## Plot the pedigree of the family with the smallest p-value.
plotPed(bir, family = "8")

}
\keyword{classes}

