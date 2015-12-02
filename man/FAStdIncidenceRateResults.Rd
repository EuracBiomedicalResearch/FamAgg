\name{FAStdIncidenceRateResults-class}
\Rdversion{1.1}
\docType{class}
\alias{$,FAStdIncidenceRateResults-method}
\alias{FAStdIncidenceRateResults-class}
\alias{FAStdIncidenceRateResults}
\alias{[,FAStdIncidenceRateResults,ANY-method}
\alias{[,FAStdIncidenceRateResults,ANY,ANY,ANY-method}
\alias{fsir,FAStdIncidenceRateResults-method}
\alias{lambda}
\alias{lambda,FAStdIncidenceRateResults-method}
\alias{plotPed,FAStdIncidenceRateResults-method}
\alias{plotRes,FAStdIncidenceRateResults-method}
\alias{result,FAStdIncidenceRateResults-method}
\alias{resultForId}
\alias{resultForId,FAStdIncidenceRateResults-method}
\alias{runSimulation,FAStdIncidenceRateResults-method}
\alias{timeInStrata}
\alias{timeInStrata,FAStdIncidenceRateResults-method}
\alias{trait<-,FAStdIncidenceRateResults-method}


\title{Familial Standardized Incidence Rate}

\description{
  The \code{FAStdIncidenceRateResults} object contains the results from a
  familial standardized incidence rate (FSIR) calculation employing in
  addition Monte Carlo simulations to assess significance levels for the
  individuals' FSIRs.

  The FSIR weights the disease status of relatives based on their degree
  of relatedness with the proband (Kerber, 1995). Formally, the FSIR is
  defined as the standardized incidence ratio (SIR) or standardized
  morality ratio in epidemiology, i.e. as the ratio between the observed
  and expected number of cases, only that both are in addition also
  weighted by the degree of relatedness (i.e. kinship value) between
  individuals in the pedigree.

}

\section{Objects from the Class}{
  \code{FAStdIncidenceRateResults} objects are created by the
  \code{\link{fsirTest}} method on a
  \code{\linkS4class{FAData}} object.
}
\usage{

\S4method{fsir}{FAStdIncidenceRateResults}(object, trait=NULL, lambda=NULL,
                                           timeInStrata=NULL,
                                           prune=TRUE, ...)

\S4method{lambda}{FAStdIncidenceRateResults}(object, ...)

\S4method{plotPed}{FAStdIncidenceRateResults}(object, id=NULL, family=NULL,
                                              filename=NULL, device="plot",
                                              only.phenotyped=FALSE, ...)

\S4method{plotRes}{FAStdIncidenceRateResults}(object, id=NULL, family=NULL,
                                              addLegend=TRUE, type="density", ...)

\S4method{result}{FAStdIncidenceRateResults}(object, method="BH")

\S4method{resultForId}{FAStdIncidenceRateResults}(object, id=NULL)

\S4method{runSimulation}{FAStdIncidenceRateResults}(object, nsim=50000,
                                                    lambda=NULL,
                                                    timeInStrata=NULL,
                                                    strata=NULL,
                                                    prune=TRUE, ...)

\S4method{timeInStrata}{FAStdIncidenceRateResults}(object)

\S4method{trait}{FAStdIncidenceRateResults}(object) <- value

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

    For \code{resultForId}: the ID of the individual from which the
    result should be returned.
  }

  \item{lambda}{
    Numeric vector with the incidence rates per stratum from the
    population. The length of this vector has to match the number of
    columns of argument \code{timeInStrata}.

    For \code{fsir}: not used.
  }

  \item{method}{
    The multiple hypothesis testing method. All methods supported by
    \code{\link{p.adjust}} are allowed.
  }

  \item{nsim}{
    Number of simulations.
  }

  \item{object}{
    The \code{FAStdIncidenceRateResults} object.
  }

  \item{only.phenotyped}{
    For \code{plotPed}: Wheter only phenotyped individuals,
    i.e. individuals with a non-NA value in column \code{affected} (the
    trait information). If \code{TRUE}, the function removes all
    non-phenotyped individuals, keeping only those that are required for
    the pedigree to be complete.
  }

  \item{prune}{
    For \code{fsir}: not used.
  }

  \item{strata}{
    For \code{runSimulation}: a numeric, character of factor
    characterizing each individual in the pedigree. The length of this
    vector and the ordering has to match the pedigree.
    This vector allows to perform stratified random sampling. See
    details on the \code{\link{PedigreeAnalysis}} help page or examples
    for more information.

  }

  \item{timeInStrata}{
    For \code{runSimulation}: a numeric matrix specifying the time at
    risk for each individual in each strata. Columns represent the
    strata, rows the individuals, each cell the time at risk for the
    individual in the respective strata. See example below. The
    \code{\link{factor2matrix}} could be useful in generating such a
    table.

    For \code{fsir}: not used.
  }

  \item{trait}{
    For \code{fsir}: not used.
  }

  \item{type}{
    For \code{plotRes}: at present only \code{"density"} is supported.
  }

  \item{value}{
    For \code{trait<-}: can be a named numeric, character or factor
    vector. The names (at least some of them) have to match the ids in
    the pedigree of the object.
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
    \item{timeInStrata}{
      Numeric matrix specifying the time at risk in each strata
      (represented by the columns) of each individual (rows). Use the
      accessor method \code{timeInStrata} to
      extract this data.
    }
    \item{lambda}{
      Numeric vector with the incidence rates per stratum from the
      population.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{object$name}{
      Access various results and data stored in the object. \code{name}
      can be \code{fsir}, \code{pvalue}, \code{lambda} or
      \code{timeInStrata} to access the FSIR, the (raw) p-value from the
      simulation analysis, the lambda or the time in strata.
    }

    \item{fsir}{
      Returns the FSIR values calculated by a call to
      the \code{runSimulation} method or
      \code{\link{fsir}} method on a
      \code{\link{FAData}} object. In contrast to that latter method,
      which directly calculates the values, this method returns the
      values from a calculation stored inside the
      \code{FAStdIncidenceRateResults} object.

      The method returns a named numeric vector with the familial
      incidence rates, the names being the ID of the individuals.
    }

    \item{plotPed}{
      Plots a pedigree for one of the affected individuals in the
      simulation results. The id of the selected affected individual
      (specified with argument \code{id}) is highlighted in red. The
      FSIR value for each individual is drawn below
      the individual's id.
      See \code{\link{plotPed}} for more details.
    }

    \item{plotRes}{
      Plots the distribution of expected FSIR
      calculated for the selected individuals from Monte Carlo
      simulations along with the actually observed FSIR.
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
      \code{"fsir"}: the familial standardized incidence rate.
      \code{"pvalue"}: the p-value for the significance of the familial
      standardized incidence rate assessed by Monte Carlo simulations.
      \code{"padj"}: the p-value adjusted for multiple hypothesis
      testing (with the method specified with argument \code{method}).

      The returned \code{data.frame} is sorted by column
      \code{"pvalue"}, its row names correspond to column
      \code{"id"}.
    }

    \item{resultForId}{
      Extracts results information for a given individual. The method
      returns a \code{list} with elements \code{"id"}, \code{"fsir"},
      \code{"pvalue"}, \code{"timeInStrata"} and \code{"lambda"} with
      the ID of the individual, the FSIR and corresponding p-value
      estimated by the simulation, the row from the \code{timeInStrata}
      matrix of the individual and the \code{lambda} (incidence rates
      from the population).
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree and trait
      information stored in the object as well as the time at risk
      provided with argument \code{timeAtRisk}.
      Returns a \code{FAStdIncidenceRateResults} object with the results
      from the simulation.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\details{
  Monte Carlo simulation and empirical p-value estimation: see details
  of \code{\link{FAIncidenceRateResults}} as the concept and calculation
  is essentially identical.

  Calling the \code{runSimulation} method on a \code{FAStdIncidenceRateResults}
  object is the same as calling the
  \code{\link{fsirTest}} on a
  \code{\linkS4class{FAData}} object. In the first case the simulation
  is performed using the trait information data stored internally in the
  object, while in the latter case the trait information have to be
  submitted to the function call.

  By providing argument \code{strata}, the stratified random sampling is
  performed. See example below and the details section in
  \code{\link{PedigreeAnalysis}} for more details.

  The FSIR can also be directly calculated, without
  simulation, using the \code{\link{fsir}} method of a
  \code{\link{FAData}} object.

  A call to the setter methods \code{trait<-} resets any simulation
  results present in the \code{sim} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.

  By default the Monte Carlo p-value estimation in the
  \code{runSimulation} method is quite memory demanding. For very large
  pedigrees the optional argument \code{lowMem=TRUE} might be passed to
  the method which results in faster and less memory demanding
  calculations. This will however disable the \code{plotRes} method on
  the resulting \code{FAStdIncidenceRateResults} as the distribution of
  FSIR from the simulation runs is no longer reported.

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

## Next we have to specify lambda and the timeInStrata matrix.
## For lambda we use information from Cancer Research UK:
## New breast cancer cases in females per year and 100000: 155.3
## New breast cancer cases in males per year and 100000:     1.1
## New prostate cancer cases in females per year and 100000: 0
## New prostate cancer cases in males per year and 100000: 134.3
lbda <- c(M=(1.1+134.4)/100000, F=155.3/100000)

## Next we need the time at risk of each individual in each strata.
## For strata we use male and females, for the time at risk we use
## column "endage":
stratMat <- factor2matrix(fad$sex)
## Next we have to multiply that with the endage, since each person
## spent that time "at risk" to get cancer.
stratMat <- stratMat * mbsub$endage

## Running the simulation
fsirs <- fsirTest(fad, trait=tcancer, traitName="cancer", lambda=lbda,
                  timeInStrata=stratMat, nsim=500)

## Showing some of the results
head(result(fsirs))

## Extract the FSIR
head(fsirs$fsir)

## We can also directly access the p-values from the simulation,
## these are however the raw, unadjusted p-values.
head(fsirs$pvalue)

## Get the time in strata for each individual
head(fsirs$timeInStrata)

## Plot the pedigree for a family with significant FSIRs.
## The numbers below the IDs of the individuals represent the actual
## FSIR values.
plotPed(fsirs, family=result(fsirs)$family[1])

## Plot also the result from the simulation run.
plotRes(fsirs, id=result(fsirs)$id[1])

## Extract the data and result for an individual
resultForId(fsirs, id="16424")

}
\keyword{classes}

