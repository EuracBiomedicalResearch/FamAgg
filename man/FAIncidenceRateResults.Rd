\name{FAIncidenceRateResults-class}
\Rdversion{1.1}
\docType{class}
\alias{$,FAIncidenceRateResults-method}
\alias{FAIncidenceRateResults-class}
\alias{FAIncidenceRateResults}
\alias{[,FAIncidenceRateResults,ANY-method}
\alias{[,FAIncidenceRateResults,ANY,ANY,ANY-method}
\alias{familialIncidenceRate,FAIncidenceRateResults-method}
\alias{plotPed,FAIncidenceRateResults-method}
\alias{plotRes,FAIncidenceRateResults-method}
\alias{result,FAIncidenceRateResults-method}
\alias{runSimulation,FAIncidenceRateResults-method}
\alias{timeAtRisk}
\alias{timeAtRisk<-}
\alias{timeAtRisk,FAIncidenceRateResults-method}
\alias{timeAtRisk<-,FAIncidenceRateResults-method}
\alias{trait<-,FAIncidenceRateResults-method}


\title{Familial Incidence Rate}

\description{
  The \code{FAIncidenceRateResults} object contains the results from a
  familial incidence rate calculation employing in addition Monte Carlo
  simulations to assess significance levels for the familial incidence
  rate of each individual.
  The familial incidence rate (FIR, also referred to as FR in Kerber
  (1995)) is an estimate for the risk per gene-time for each individuals
  for a certain disease given the disease experience in the cohort. The
  measure considers the kinship of each individual with any affected
  individual in the pedigree and the time at risk for each.

  Note that in contrast to e.g. \code{\link{FAKinSumResults}} a familial
  incidence rate and corresponding p-value are calculated and available
  for all individuals in the pedigree, not only for affected individuals.
}

\section{Objects from the Class}{
  \code{FAIncidenceRateResults} objects are created by the
  \code{\link{familialIncidenceRateTest}} method on a
  \code{\linkS4class{FAData}} object.
}
\usage{

\S4method{familialIncidenceRate}{FAIncidenceRateResults}(object, trait=NULL,
                                                         timeAtRisk=NULL,
                                                         ...)

\S4method{plotPed}{FAIncidenceRateResults}(object, id=NULL, family=NULL,
                                           filename=NULL, device="plot",
                                           only.phenotyped=FALSE, ...)

\S4method{plotRes}{FAIncidenceRateResults}(object, id=NULL, family=NULL,
                                           addLegend=TRUE, type="density", ...)

\S4method{result}{FAIncidenceRateResults}(object, method="BH")

\S4method{runSimulation}{FAIncidenceRateResults}(object, nsim=50000,
                                                 timeAtRisk=NULL,
                                                 strata=NULL, ...)

\S4method{timeAtRisk}{FAIncidenceRateResults}(object)

\S4method{timeAtRisk}{FAIncidenceRateResults}(object) <- value

\S4method{trait}{FAIncidenceRateResults}(object) <- value

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
    For \code{plotPed}: the family for which the pedigree should be
    plotted.
    For \code{plotRes}: not supported.
  }

  \item{filename}{
    For \code{plotPed}: the file name to which the pedigree plot should
    be exported. See \code{\link{plotPed}} for more details.
  }

  \item{id}{
    For \code{plotPed} and \code{plotRes}: the id of the indiviual
    for which the pedigree or the simulation result should be plotted.
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
    The \code{FAIncidenceRateResults} object.
  }

  \item{only.phenotyped}{
    For \code{plotPed}: Wheter only phenotyped individuals,
    i.e. individuals with a non-NA value in column \code{affected} (the
    trait information). If \code{TRUE}, the function removes all
    non-phenotyped individuals, keeping only those that are required for
    the pedigree to be complete.
  }

  %% \item{rm.singletons}{
  %%   For \code{familialIncidenceRate}: not used, for
  %%   \code{familialIncidenceRateTest}: remove singletons from the data
  %%   set prior to the analysis. This does however not have an impact on
  %%   the results, as the FIR for singletons will anyway be \code{NA}
  %%   because they do not share kinship with any other individual in the
  %%   pedigree.
  %% }

  \item{strata}{
    For \code{runSimulation}: a numeric, character of factor
    characterizing each individual in the pedigree. The length of this
    vector and the ordering has to match the pedigree.
    This vector allows to perform stratified random sampling. See
    details on the \code{\link{PedigreeAnalysis}} help page or examples
    for more information.

  }

  \item{timeAtRisk}{
    For \code{runSimulation}: a numeric vector specifying the time at
    risk for each individual in the given trait. See also
    \code{\link{estimateTimeAtRisk}}, an utility function to estimate
    time at risk.
    For \code{familialIncidenceRate}: not used.

  }

  \item{trait}{
    For \code{familialIncidenceRate}: not used.
  }

  \item{type}{
    For \code{plotRes}: at present only \code{"density"} is supported.
  }

  \item{value}{
    For \code{trait<-}: can be a named numeric, character or factor
    vector. The names (at least some of them) have to match the ids in
    the pedigree of the object.
    For \code{timeAtRisk}: a numeric vector with the time at risk for
    each individual.
  }

  \item{...}{
    For \code{plotPed}: additional arguments to be submitted to the
    internal \code{buildPed} call and to \code{\link{plotPed}}.

    For \code{runSimulation}: additional arguments \code{prune} and
    \code{lowMem}. See below for details.
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
      directly, use the \code{result} method to extract result
      information.
    }
    \item{timeAtRisk}{
      Numeric vector with the time at risk for each individual. Use the
      accessor method \code{timeAtRisk} or use \code{object$tar} to
      extract this data.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{object$name}{
      Access the familial incidence rate using \code{object$fir}, the
      (raw) p-value from the simulation using \code{object$pvalue} and the
      time at risk for each individual using \code{object$tar} or
      \code{object$timeAtRisk} with \code{object} being the
      \code{FAIncidenceRateResults} object.
    }

    \item{familialIncidenceRate}{
      Returns the familial incidence rate values calculated by a call to
      the \code{runSimulation} method or
      \code{\link{familialIncidenceRate}} method on a
      \code{\link{FAData}} object. In contrast to that latter method,
      which directly calculates the values, this method returns the
      values from a calculation stored inside the
      \code{FAIncidenceRateResults} object.

      The method returns a named numeric vector with the familial
      incidence rates for all individuals in the pedigree, the names
      being the ID of the individuals. Singletons as well as individuals
      that, after removing not phenotyped individuals or individuals
      without time at risk, do not share kinship with any other
      individual in the pedigree have a value of \code{NA}.
    }

    \item{plotPed}{
      Plots a pedigree for one of the affected individuals in the
      simulation results. The id of the selected affected individual
      (specified with argument \code{id}) is highlighted in red. The
      familial incidence rate value for each individual is drawn below
      the individual's id.
      See \code{\link{plotPed}} for more details.
    }

    \item{plotRes}{
      Plots the distribution of expected familial incidence rates
      calculated for the selected individuals from Monte Carlo
      simulations along with the actually observed familial incidence
      rate.
    }

    \item{result}{
      Returns the result from the simulation as a \code{data.frame} with
      columns:

      \code{"trait_name"}: the name of the trait.

      \code{"total_phenotyped"}: total number of phenotyped individuals in
      the trait.

      \code{"total_affected"}: total number of affected individuals in the
      trait.

      \code{"total_tested"}: the number of individuals in the pedigree
      considered for the simulation. This corresponds to all individuals
      with valid, non-NA, values in trait, timeAtRisk and eventually
      strata.

      \code{"id"}: the id of the individual.

      \code{"family"}: the family id.

      \code{"fir"}: the familial incidence rate. Note that this will be

      \code{NA} for all non-phenotyped individuals and singletons in the pedigree as well as for
      individuals that do not share kinship with at least one other
      phenotyped individual with valid time at risk (or valid value in
      parameter \code{strata}).

      \code{"pvalue"}: the p-value for the significance of the familial
      incidence rate assessed by Monte Carlo simulations.

      \code{"padj"}: the p-value adjusted for multiple hypothesis
      testing (with the method specified with argument \code{method}).

      The returned \code{data.frame} is sorted by column
      \code{"pvalue"}, its row names correspond to column
      \code{"id"}.
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree and trait
      information stored in the object as well as the time at risk
      provided with argument \code{timeAtRisk}.
      Returns a \code{FAIncidenceRateResults} object with the results
      from the simulation.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\details{
  Monte Carlo simulation and empirical p-value estimation: the
  background distribution to calculate the p-value for a familial
  incidence rate (FIR) is determined by randomly sampling N affected
  individuals (N being the number of affected) and calculating the
  \emph{expected} FIR for all individuals in each simulation
  iteration. The p-value for an individual represents thus the number of
  times an expected FIR for that individual from the simulation was
  found to be larger than or equal to the observed FIR divided by the
  number of iterations.

  Calling the \code{runSimulation} method on a \code{FAIncidenceRateResults}
  object is the same as calling the
  \code{\link{familialIncidenceRateTest}} on a
  \code{\linkS4class{FAData}} object. In the first case the simulation
  is performed using the trait information data stored internally in the
  object, while in the latter case the trait information have to be
  submitted to the function call.

  By providing argument \code{strata}, the stratified random sampling is
  performed. See example below and the details section in
  \code{\link{PedigreeAnalysis}} for more details.

  The familial incidence rate can also be directly calculated, without
  simulation, using the \code{\link{familialIncidenceRate}} method of a
  \code{\link{FAData}} object.

  A call to the setter methods \code{trait<-} resets any simulation
  results present in the \code{sim} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.

  By default the Monte Carlo p-value estimation in the
  \code{runSimulation} method is quite memory demanding. For very large
  pedigrees the optional argument \code{lowMem=TRUE} might be passed to
  the method which results in faster and less memory demanding
  calculations. This will however disable the \code{plotRes} method on
  the resulting \code{FAIncidenceRateResults} as the distribution of
  familial incidence rates from the simulation runs is no longer
  reported.

  Note: the FIR for singletons and individuals that do not share kinship
  with at least one other phenotyped individual that has also a valid
  value in argument \code{timeAtRisk} (and eventually \code{strata})
  will be \code{NA}.
}
\note{
  Subsetting (using the \code{[} operator) is not supported.
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer
}
\references{
  Kerber, R.A. (1995) Method for calculating risk associated with family
  history of a disease. \emph{Genet Epidemiol}, pp 291--301.
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
  \code{\link{plotPed}},
  \code{\link{estimateTimeAtRisk}}
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
                    minnbreast$famid == 6 | minnbreast$famid == 7 |
                    minnbreast$famid == 411, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- c("family", "id", "father", "mother", "sex")

## Generate the FAData.
fad <- FAData(pedigree=PedDf)

## Specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## Spefify the "time at risk"; we are using column "endage"
tar <- mbsub$endage

## Perform the simulation test:
far <- familialIncidenceRateTest(fad, trait=tcancer, traitName="cancer",
                                 timeAtRisk=tar, nsim=1000)
head(result(far))

## We can easily extract the actual FIR values:
head(far$fir)
## Or
head(familialIncidenceRate(far))

## Access the p-value directly.
head(far$pvalue)

## Access the time at risk
head(timeAtRisk(far))
head(far$tar)
head(far$timeAtRisk)

## Plot the pedigree for a family with significant FIRs.
## The numbers below the IDs of the individuals represent the actual
## FIR values.
plotPed(far, family=result(far)$family[1])

## Plot also the result from the simulation run.
plotRes(far, id=result(far)$id[1])

}
\keyword{classes}

