\name{FAKinSumResults-class}
\Rdversion{1.1}
\docType{class}
\alias{FAKinSumResults-class}
\alias{FAKinSumResults}
\alias{[,FAKinSumResults,ANY-method}
\alias{[,FAKinSumResults,ANY,ANY,ANY-method}
\alias{plotPed,FAKinSumResults-method}
\alias{plotRes,FAKinSumResults-method}
\alias{result,FAKinSumResults-method}
\alias{runSimulation,FAKinSumResults-method}
\alias{trait<-,FAKinSumResults-method}


\title{Kinship sum test}

\description{
  The \code{FAKinSumResults} object contains the results from a kinship
  cluster test which evaluates familial aggregation based on the sum of kinship
  values between affected cases. This test highlights individuals that
  exhibit a higher than chance relationship to other affected
  individuals, therefore highlighting individuals within families
  aggregating the phenotype. To achieve this, for each affected
  individual the sum of kinship values to all other affected cases is
  computed. In a Monte Carlo simulation this is repeated with the same
  number of cases and the resulting background distribution is used to
  compute p-values for the kinship sums obtained from the observed
  cases.
}
\section{Objects from the Class}{
  \code{FAKinSumResults} objects are created by the
  \code{\link{kinshipSumTest}} method on a
  \code{\linkS4class{FAData}} object.
}
\usage{

\S4method{plotPed}{FAKinSumResults}(object, id=NULL, family=NULL,
                                     filename=NULL, device="plot",
                                     only.phenotyped=FALSE, ...)

\S4method{plotRes}{FAKinSumResults}(object, id=NULL, family=NULL,
                                     addLegend=TRUE, type="density", ...)

\S4method{result}{FAKinSumResults}(object, method="BH")

\S4method{runSimulation}{FAKinSumResults}(object, nsim=50000,
strata=NULL, ...)

\S4method{trait}{FAKinSumResults}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{addLegend}{
    For \code{plotRes}: if a legend should be added to the plot.
  }

  \item{device}{
    For \code{plotPed}: see \code{\link{plotPed}} for more details.
  }

  \item{family}{
    For \code{plotPed}: not supported.
  }

  \item{filename}{
    For \code{plotPed}: the file name to which the pedigree plot should
    be exported. See \code{\link{plotPed}} for more details.
  }

  \item{id}{
    For \code{plotPed} and \code{plotRes}: the id of the indiviual
    (i.e. affected individual in the result \code{data.frame}) for which
    the pedigree or the simulation result should be plotted.
    Note: \code{id} can be a \code{numeric} or a
    \code{character}. Numeric ids will be internally converted to
    \code{character}.
  }

  \item{method}{
    The multiple hypothesis testing method. All methods supported by
    \code{\link{p.adjust}} are allowed.
  }

  \item{nsim}{
    Number of simulations.
  }

  \item{object}{
    The \code{FAKinSumResults} object.
  }

  \item{only.phenotyped}{
    For \code{plotPed}: Wheter only phenotyped individuals,
    i.e. individuals with a non-NA value in column \code{affected} (the
    trait information). If \code{TRUE}, the function removes all
    non-phenotyped individuals, keeping only those that are required for
    the pedigree to be complete.
  }

  \item{strata}{
    For \code{runSimulation}: a numeric, character of factor
    characterizing each individual in the pedigree. The length of this
    vector and the ordering has to match the pedigree.
    This vector allows to perform stratified random sampling. See
    details on the \code{\link{PedigreeAnalysis}} help page or examples
    for more information.

  }

  \item{type}{
    For \code{plotRes}: either \code{"density"} (the default) or
    \code{"hist"} specifying whether the distribution of expected values
    from the simulation should be visualized as a density plot or
    histogram.
  }

  \item{value}{
    For \code{trait<-}: can be a named numeric, character or factor
    vector. The names (at least some of them) have to match the ids in
    the pedigree of the object.
  }

  \item{...}{
    For \code{plotPed}: additional arguments to be submitted to the
    internal \code{buildPed} call and to \code{\link{plotPed}}.
  }

}
\section{Extends}{
  Class \code{\linkS4class{FAData}} directly.
}
\section{Slots}{
  \describe{
    \item{nsim}{
      Number of simulations.
    }
    \item{sim}{
      The result of the simulation. This slot should not be accessed
      directly, use the \code{result} method to extract result information.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{plotPed}{
      Plots a pedigree for one of the affected individuals in the
      simulation results. The id of the selected affected individual
      (specified with argument \code{id}) is highlighted in red.
      See \code{\link{plotPed}} for more details.
    }

    \item{plotRes}{
      Plots the distribution of kinship sums between random sets of
      samples individuals from the Monte Carlo simulation along with the
      actually observed kinship sum for the affected individual
      specified with parameter \code{id}. For \code{id} only affected
      individuals for which the analysis has been performed are
      allowed. The ids of these individuals are listed in column
      \code{"affected_id"} of the \code{data.frame} returned by
      \code{result}.

    }

    \item{result}{
      Returns the result from the simulation as a \code{data.frame} with
      columns:
      \code{"trait_name"}: the name of the trait.
      \code{"total_phenotyped"}: total number of phenotyped individuals in
      the trait.
      \code{"total_affected"}: total number of affected individuals in the
      trait.
      \code{"affected_id"}: the id of the affected individual for whom
      the test has been performed.
      \code{"family"}: the family id of the affected.
      \code{"kinship_sum"}: the sum of kinship values.
      \code{"freq"}: the expected frequency of the kinship sum from the
      simulation.
      \code{"pvalue"}: the p-value for the significance of the kinship
      sum.
      \code{"padj"}: the p-value adjusted for multiple hypothesis
      testing (with the method specified with argument \code{method}).

      The returned \code{data.frame} is sorted by column
      \code{"pvalue"}, its row names correspond to column
      \code{"affected_id"}.
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree and trait
      information stored in the object. Returns a
      \code{FAKinSumResults} object with the results of the simulation.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\details{
  Calling the \code{runSimulation} method on a \code{FAKinSumResults}
  object is the same as calling the \code{\link{kinshipSumTest}} on a
  \code{\linkS4class{FAData}} object. In the first case the
  simulation is performed using the trait information
  data stored internally in the object, while in the latter case the
  trait information have to be submitted to the function
  call.

  A call to the setter methods \code{trait<-} resets any simulation
  results present in the \code{sim} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.

  The expected frequency (column \code{"freq"}) in the result
  \code{data.frame} is \code{NA} if the corresponding sum of kinship
  coefficients reported in column \code{"kinship_sum"} was never sampled
  in the simulation. Still, a p-value can be reported.

  The \code{plotPed} function does not support to draw pedigrees for
  individuals for which no simulation test has been performed. To draw a
  pedigree for any individual (with or without trait information, being
  affected or not in the trait) refer to the \code{plotPed} method for
  \code{\linkS4class{FAData}} objects.
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
##  Perform the simulation analysis
##
## Load the test data.
data(minnbreast)

## Subset to some families and generate a pedigree data.frame.
mbsub <- minnbreast[minnbreast$famid == 4 | minnbreast$famid == 5 |
                  minnbreast$famid == 6 | minnbreast$famid == 7, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- c("family", "id", "father", "mother", "sex")

## Generate the FAData.
fad <- FAData(pedigree=PedDf)

## Specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## Perform the test:
far <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                      nsim=1000)
head(result(far))

## Plot the pedigree for one of affected individuals. The id of the affected
## individual is highlighted in red.
plotPed(far, id=result(far)$affected_id[1])

## Replace the trait, this will delete all simulation results
## on the existing FAKinSumResults object
tpreg <- mbsub$everpreg
names(tpreg) <- mbsub$id
trait(far) <- tpreg

## The analysis can be performed using the runSimulation method.
far <- runSimulation(far, nsim=1000)
head(result(far))

## Plot the pedigree of one of the affected; note that "affected" now
## indicates whether the individual was ever pregnant.
plotPed(far, id="9")

## Plot also the result from the simulation run.
plotRes(far, id="9")

}
\keyword{classes}

