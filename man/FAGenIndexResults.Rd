\name{FAGenIndexResults-class}
\Rdversion{1.1}
\docType{class}
\alias{FAGenIndexResults-class}
\alias{FAGenIndexResults}
\alias{[,FAGenIndexResults,ANY-method}
\alias{[,FAGenIndexResults,ANY,ANY,ANY-method}
\alias{plotPed,FAGenIndexResults-method}
\alias{plotRes}
\alias{plotRes,FAGenIndexResults-method}
\alias{result,FAGenIndexResults-method}
\alias{runSimulation,FAGenIndexResults-method}
\alias{trait<-,FAGenIndexResults-method}


\title{Genealogical Index}

\description{
  The genealogical index (Hill, 1980) is a method to identify
  familial clustering of diseases or other traits. For a given trait,
  the method computes the mean kinship between affected in the
  whole pedigree along with mean kinships of randomly drawn sets of
  individuals. The distribution of average kinship values among the
  control sets is used to estimate the probability that the observed
  level of kinship among the cases is due to chance.
}
\details{
  This implementation differs from the original method from Hill as it
  allows, in addition to perform per family analyses, to use also
  stratified sampling and allows a more flexible definition of the set
  of matched control individuals.
  The \code{controlSetMethod} parameter allows to specify a method to
  define the matched control set (e.g. matched by sex or matched by any
  externally provided vector).

  Stratified sampling allows to even further fine tune the selection of
  matched controls. Assuming that in a pedigree the group of affected
  consists of 5 females and 3 male individuals, passing the sex of all
  individuals to the function (e.g. \code{strata=fad$sex}, with
  \code{fad} being the \code{\linkS4class{FAData}} object containing the
  pedigree to be analyzed) results in random sets with the same
  proportion of male/female individuals (i.e. 5 females, 3 males).

  Note that, if \code{strata} is specified, all individuals with a
  missing value in \code{strata} (also affected individuals) are
  excluded from the analysis.

  By default, the genealogical index is calculated on the whole
  pedigree, but it is also possible to evaluate within-family clustering
  of cases by specifying \code{perFamilyTest=TRUE}. In that case, it is
  also possible to use the \code{\link{getGenerationMatched}} and
  \code{\link{getGenerationSexMatched}} functions to define the set of
  matched controls from which random samples will be taken.

  A call to the setter methods \code{trait<-} resets any simulation
  results present in the \code{sim} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.


}
\section{Objects from the Class}{
  \code{FAGenIndexResults} objects are created calling the
  \code{\link{genealogicalIndexTest}} method on a
  \code{\linkS4class{FAData}} object.
}
\usage{

\S4method{plotPed}{FAGenIndexResults}(object, id=NULL, family=NULL,
                                     filename=NULL, device="plot", ...)

\S4method{plotRes}{FAGenIndexResults}(object, id=NULL, family=NULL,
                                     addLegend=TRUE, type="density", ...)

\S4method{result}{FAGenIndexResults}(object, method="BH")

\S4method{runSimulation}{FAGenIndexResults}(object, nsim=50000,
                                           perFamilyTest=FALSE,
                                           controlSetMethod="getAll",
                                           prune=TRUE, strata=NULL, ...)

\S4method{trait}{FAGenIndexResults}(object) <- value

}
\arguments{

  (in alphabetic order)

  \item{addLegend}{
    For \code{plotRes}: if a legend should be added to the plot.
  }

  \item{controlSetMethod}{
    For \code{runSimulation}: the method (i.e. name of the function)
    that should be used to define the set of (eventually matched)
    control individuals from which the random samples are
    taken. Supported functions are \code{\link{getAll}},
    \code{\link{getSexMatched}} and
    \code{\link{getExternalMatched}}. For \code{perFamilyTest=TRUE} also
    \code{\link{getGenerationMatched}} and
    \code{\link{getGenerationSexMatched}} are supported. Note: for
    \code{getExternalMatched}, a numeric, character or factor vector to
    be used for the matching has to be submitted to \code{runSimulation}
    as additional argument \code{match.using}.
  }

  \item{device}{
    For \code{plotPed}: see \code{\link{plotPed}} for more details.
  }

  \item{family}{
    For \code{plotPed}: the family for which the pedigree should be
    plotted.
    For \code{plotRes}: the family for which the genealogical index
    analysis simulation results should be shown. Only supported if
    \code{perFamilyTest=TRUE}.
  }

  \item{filename}{
    For \code{plotPed}: the file name to which the pedigree plot should
    be exported. See \code{\link{plotPed}} for more details.
  }

  \item{id}{
    For \code{plotPed}: the id of an indiviual from a family for which
    the pedigree should be plotted.
    For \code{plotRes}: the id of an individual from a family for which
    the genealogical index analysis simulation results should be
    shown. Only supported if \code{perFamilyTest=TRUE}.
  }

  \item{method}{
    The multiple hypothesis testing method. All methods supported by
    \code{\link{p.adjust}} are allowed.
  }

  \item{nsim}{
    Number of simulations.
  }

  \item{perFamilyTest}{
    For \code{runSimulation}: whether the test should be performed
    on the whole pedigree (default) or separately within each
    family. In the latter case the test evaluates the presence of
    clustered affected individuals within each family.
  }

  \item{prune}{
    For \code{runSimulation}: whether unconnected individuals in the
    pedigree should be removed.
  }

  \item{object}{
    The \code{FAGenIndexResults} object.
  }

  \item{strata}{
    For \code{runSimulation}: a numeric, character of factor
    characterizing each individual in the pedigree. The length of this
    vector and the ordering has to match the pedigree.
    This vector allows to perform stratified random sampling. See
    details or examples for more information.
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

    For \code{runSimulation}: additional arguments passed to the choosen
    \code{controlSetMethod} function (e.g. \code{match.using} for
    \code{\link{getExternalMatched}}).
  }

}
\section{Extends}{
  Class \code{\linkS4class{FAData}} directly.
}
\section{Slots}{
  \describe{
    \item{controlSetMethod}{
      A character specifying the name of the method used to define the
      set of control individuals from which random samples were taken.
    }
    \item{nsim}{
      Number of simulations.
    }
    \item{perFamilyTest}{
      Logical indicating whether a per-family test was performed.
    }
    \item{sim}{
      The result of the simulation. This slot should not be accessed
      directly, use the \code{result} method to extract result information.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \code{plotPed}{
      Plots a pedigree for one of the affected individuals in the
      simulation results. The id of the selected affected individual
      (specified with argument \code{id}) is highlighted in red.
      See \code{\link{plotPed}} for more details.
    }

    \code{plotRes}{
      Plots the results from a genealogical index simulation
      analysis. The distribution of the mean kinship values of the
      randomly drawn controls are displayed as a grey density plot, the
      observed mean kinship value of all affected as a blue vertical
      line.
    }

    \item{result}{
      Returns the result from the simulation as a \code{data.frame} with
      columns:
      \code{"trait_name"}: the name of the trait.
      \code{"total_phenotyped"}: total number of individuals in the
      pedigree phenotyped in the analyzed trait.
      \code{"total_affected"}: total number ofindividuals in the
      pedigree that are affected in the analyzed trait (i.e. number of
      cases).
      \code{"entity_id"}: the id for the analyzed entity, being either
      the whole pedigree (in which case the id will be \code{"1"}) or
      the id of the family (if \code{perFamilyTest=TRUE}).
      \code{"entity_ctrls"}: the number of (matched) control individuals
      from which the random samples were drawn.
      \code{"entity_affected"}: the number of affected individuals in
      the entity. This number can differ from the number of affected, if
      \code{strata} was specified and some of the affected have a
      missing value in \code{strata}.
      \code{"mean_kinship"}: the mean kinship value between all affected
      in the entity.
      \code{"pvalue"}: the p-value for the significance of the mean
      kinship.
      \code{"padj"}: the p-value adjusted for multiple hypothesis
      testing (with the method specified with argument \code{method}).

      The returned \code{data.frame} is sorted by column
      \code{"pvalue"}, its rownames correspond to column
      \code{"entity_id"}.
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree and trait
      information stored in the object. Returns a
      \code{FAGenIndexResults} object with the results of the simulation.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\note{
  Subsetting (using the \code{[} operator) is not supported.
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\references{
  Hill, J. R. (1980) A survey of cancer sites by kinship in the Utah Mormon
  population. In Cairns J, Lyon JL, Skolnick M (eds):
  \emph{Cancer Incidence in Defined Populations. Banbury Report 4.}
  Cold Spring Harbor, NY: Cold Spring Harbor Laboratory Press, pp 299--318.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{FAData}},
  \code{\link{trait}},
  \code{\link{probabilityTest}},
  \code{\link{kinshipGroupTest}},
  \code{\link{kinshipSumTest}},
  \code{\link{familialIncidenceRateTest}},
  \code{\link{fsirTest}},
  \code{\link{plotPed}}
}
\examples{
##########################
##
##  Perform the simulation analysis
##
## Load the Minnesota Breast Cancer data set.
data(minnbreast)

## Subset to some families and generate a pedigree data.frame
mbsub <- minnbreast[minnbreast$famid == 4 | minnbreast$famid == 14 |
                  minnbreast$famid == 6 | minnbreast$famid == 8, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- c("family", "id", "father", "mother", "sex")

## Generate the FAData.
fad <- FAData(pedigree=PedDf)

## Specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## Perform the test with default settings, i.e. use all individuals
## in the pedigree as control set from which random samples are drawn
## and perform the analysis on the whole pedigree.
gi <- genealogicalIndexTest(fad, trait=tcancer, traitName="cancer",
                            nsim=1000,)
## Just show some information
gi

## Show the results
result(gi)

## Plot the observed mean kinship and the distribution of the mean kinship of
## random samples.
plotRes(gi)

## Plot the pedigree for one of the families. All individuals
## used as matched control set are highlighted in red.
plotPed(gi, family="8")

## Repeat the analysis using the sex as strata. This will result in stratified
## random sampling with the number of female and male individuals selected in
## each permutation corresponding to the numbers below
table(gi$sex[affectedIndividuals(gi)])
giStrata <- runSimulation(gi, nsim=1000, strata=gi$sex)
result(giStrata)


## Alternatively, we can use "getSexMatched" as the function to define the set
## of control individuals. Just, in the present case both male and females
## individuals will be selected since also there are male and female individuals
## among the affected cases.
giPerFam <- runSimulation(gi, nsim=1000, controlSetMethod="getSexMatched",
                          perFamilyTest=TRUE)
result(giPerFam)

## For those families in which there are only female cases, random samples
## were drawn among only female individuals (within the same family). These
## are highlighted in red in the pedigree plot:
plotPed(giPerFam, family="14", cex=0.5)

## Plot the simulation result for this family:
plotRes(giPerFam, family="14")

}
\keyword{classes}

