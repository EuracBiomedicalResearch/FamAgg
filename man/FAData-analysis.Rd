\name{PedigreeAnalysis}
\Rdversion{1.1}
\alias{PedigreeAnalysis}
\alias{estimateTimeAtRisk}
\alias{factor2matrix}
\alias{familialIncidenceRate}
\alias{familialIncidenceRate,FAData-method}
\alias{familialIncidenceRateTest}
\alias{familialIncidenceRateTest,FAData-method}
\alias{fsir}
\alias{fsir,FAData-method}
\alias{fsirTest}
\alias{fsirTest,FAData-method}
\alias{genealogicalIndexTest}
\alias{genealogicalIndexTest,FAData-method}
\alias{kinshipGroupTest}
\alias{kinshipGroupTest,FAData-method}
\alias{kinshipSumTest}
\alias{kinshipSumTest,FAData-method}
\alias{probabilityTest}
\alias{probabilityTest,FAData-method}
\alias{sliceAge}
\alias{binomialTest}

\title{Pedigree analysis and familial aggregation methods}

\description{
  Various functions to perform pedigree analyses and to investigate
  familial clustering of e.g. cancer cases.
}
\usage{

binomialTest(object, trait, traitName, global = FALSE, prob = NULL,
    alternative = c("greater", "less", "two.sided"))

estimateTimeAtRisk(startDate=NULL, startDateFormat="\%Y-\%m-\%d",
                   endDate=NULL, endDateFormat="\%Y-\%m-\%d",
                   incidenceDate=NULL, incidenceDateFormat="\%Y-\%m-\%d",
                   deathDate=NULL, deathDateFormat="\%Y-\%m-\%d",
                   allowNegative=FALSE, affected=NULL,
                   incidenceSubtract=0.5)

factor2matrix(x)

\S4method{familialIncidenceRate}{FAData}(object, trait=NULL,
                                         timeAtRisk=NULL)

\S4method{familialIncidenceRateTest}{FAData}(object, trait=NULL,
                                            nsim=50000, traitName=NULL,
                                            timeAtRisk=NULL,
                                            strata=NULL, ...)

\S4method{fsir}{FAData}(object, trait=NULL, lambda=NULL, timeInStrata=NULL)

\S4method{fsirTest}{FAData}(object, trait=NULL, nsim=50000, traitName=NULL,
                            lambda=NULL, timeInStrata=NULL,
                            strata=NULL, ...)

\S4method{genealogicalIndexTest}{FAData}(object, trait, nsim=50000,
                                         traitName, perFamilyTest=FALSE,
                                         controlSetMethod="getAll",
                                         rm.singletons=TRUE, strata=NULL, ...)

\S4method{kinshipGroupTest}{FAData}(object, trait, nsim=50000,
                                    traitName, strata=NULL, ...)

\S4method{kinshipSumTest}{FAData}(object, trait, nsim=50000,
                                  traitName, strata=NULL, ...)

\S4method{probabilityTest}{FAData}(object, trait, cliques,
                                   nsim=50000, traitName,
                                   ...)

sliceAge(x, slices=c(0, 40, Inf))
%connectedSubgraph(graph, nodes, mode="all", all.nodes=TRUE, ifnotfound)


}
\arguments{

  (in alphabetic order)

  \item{affected}{
    For \code{estimateTimeAtRisk}: optional parameter specifying which
    of the individuals are affected. This is useful if only \code{endDate} is
    specified, but not the \code{incidenceDate}. See method description
    for further details.
  }

  \item{allowNegative}{
    For \code{estimateTimeAtRisk}: if \code{FALSE} any negative time
    periods are set to \code{0}.
  }

  \item{alternative}{
    For \code{binomialTest}: the alternative hypothesis. See
    \code{\link{binom.test}} for more details. Defaults to
    \code{"greater"}, i.e. tests whether in a family a larger number of
    affected is present than expected by chance (given a global
    probability).
  }
  
  \item{cliques}{
    A named numeric or characted vector or factor with the names
    corresponding to ids of the individuals in the pedigree. The ids
    will be internally matched and sub-set to the ids available in the
    pedigree.
  }

  \item{controlSetMethod}{
    For \code{genealogicalIndexTest}: the method (i.e. name of the function)
    that should be used to define the set of (eventually matched)
    control individuals from which the random samples are
    taken. Supported functions are \code{\link{getAll}},
    \code{\link{getSexMatched}} and
    \code{\link{getExternalMatched}}. For \code{perFamilyTest=TRUE} also
    \code{\link{getGenerationMatched}} and
    \code{\link{getGenerationSexMatched}} are supported.
    Note: for \code{getExternalMatched}, a numeric, character or factor vector to
    be used for the matching has to be submitted
    as additional argument \code{match.using}.
  }

  \item{deathDate}{
    For \code{estimateTimeAtRisk}: the date of death.
  }

  \item{deathDateFormat}{
    For \code{estimateTimeAtRisk}: the format in which the dates are
    submitted. See \code{\link{as.Date}} for more information.
  }

  \item{endDate}{
    For \code{estimateTimeAtRisk}: the end date, which can be the end
    date for the study or, if \code{deathDate} and \code{incidenceDate}
    are not specified, the earliest time point of: date of incidence,
    death or end of study.
  }

  \item{endDateFormat}{
    For \code{estimateTimeAtRisk}: the format in which the dates are
    submitted. See \code{\link{as.Date}} for more information.
  }

  \item{global}{
    For \code{binomialTest}: whether the binomial test should be applied
    to the whole pedigree, or family-wise (default). If
    \code{global = TRUE} the population probability has to be provided
    with parameter \code{prob}.
  }
  
  \item{incidenceDate}{
    For \code{estimateTimeAtRisk}: the date of the incidence for an
    individual, i.e. the date when the status was changed from
    un-affected to affected in the to be analyzed trait.
  }

  \item{incidenceDateFormat}{
    For \code{estimateTimeAtRisk}: the format in which the dates are
    submitted. See \code{\link{as.Date}} for more information.
  }

  \item{incidenceSubtract}{
    For \code{estimateTimeAtRisk}: the amount of time (of the time unit
    of the time at risk) that should be
    subtracted from the calculated time at risk for affected
    individuals. See method description below for more details.
  }

  \item{lambda}{
    Numeric vector with the incidence rates per stratum from the
    population. The length of this vector has to match the number of
    columns of argument \code{timeInStrata}.
  }

  \item{nsim}{
    The number of simulations.
  }

  \item{object}{
    The \code{FAData} object.
  }

  \item{perFamilyTest}{
    For \code{genealogicalIndexTest}: whether the test should be performed
    on the whole pedigree (default) or separately within each
    family. In the latter case the test evaluates the presence of
    clustered affected individuals within each family.
  }

  \item{prob}{
    For \code{binomialTest}: the hypothesized probability of success
    (being affected) from/for the whole population.
  }
  
  \item{rm.singletons}{
    For \code{genealogicalIndexTest}:
    whether unconnected individuals in the pedigree (singletons) should
    be removed from the pedigree prior to the analysis.
  }

  \item{slices}{
    For \code{sliceAge}: a numeric vector defining the
    age-slices. Similar to argument \code{vec} for
    \code{\link{findInterval}}. Defines the minimum and maximum age for
    the age slices, i.e. first number corresponds to the lower boundary
    of the first age slice, the second number to the upper boundary of
    the first and lower boundary of the second age slice and so on.
  }

  \item{startDate}{
    For \code{estimateTimeAtRisk}: the date of the start of the
    study. Can also be the birth date.
  }

  \item{startDateFormat}{
    For \code{estimateTimeAtRisk}: the format in which the dates are
    submitted. See \code{\link{as.Date}} for more information.
  }

  \item{strata}{
    For \code{genealogicalIndexTest}, \code{kinshipGroupTest} and
    \code{kinshipSumTest}: a numeric, character or factor
    characterizing each individual in the pedigree. The length of this
    vector and the ordering has to match the pedigree.
    This vector allows to perform stratified random sampling. See
    details for more information.
  }

  \item{timeAtRisk}{
    A numeric vector specifying the \emph{time at risk} for each
    individual. The definition for this variable is taken from Kerber
    (1995). See description of the method below for more information.
    \code{timeAtRisk} has to have the same number of elements than there are
    individuals in the pedigree and it is assumed that the ordering of
    the vector matches the order of the individuals in the pedigree.
  }

  \item{timeInStrata}{
    For \code{fsir} and \code{fsirTest}: a numeric matrix specifying the time at
    risk for each individual in each strata. Columns represent the
    strata, rows the individuals, each cell the time at risk for the
    individual in the respective strata.
  }

  \item{trait}{
    A named numeric vector (values \code{0}, \code{1} and \code{NA}) or
    logical vector (values \code{FALSE}, \code{TRUE} and \code{NA}) with
    the names matching the ids of the individuals in the pedigree. The
    method internally matches and re-orders the trait vector to match
    the ordering of the ids in the pedigree.

    If trait is not specified, the trait information stored within the
    \code{FAData} object is used.
  }

  \item{traitName}{
    The name of the trait (optional).
  }

  \item{x}{
    For \code{sliceAge}: a numeric vector representing the age of
    individuals.
    For \code{factor2matrix}: a factor that should be converted into a
    matrix.
  }

  \item{...}{
    For \code{genealogicalIndexTest}: additional arguments passed to the choosen
    \code{controlSetMethod} function (e.g. \code{match.using} for
    \code{\link{getExternalMatched}}).

    For \code{familialIncidenceRateTest}: use \code{lowMem=TRUE} for
  very large pedigrees. This will use a faster and less memory demanding
  p-value estimation.
  }

}
\details{
  Stratified sampling: some of the familial aggregation methods allow to
  use stratified sampling for the Monte Carlo simulations. In stratified
  sampling, the same number of random samples will be selected within each
  class/stratum then there are among the affected. As example, if 5
  female and 2 male individuals are affected in the analysed trait and
  sex stratified sampling is performed, in each permuatation the same
  number of random samples in each group (i.e. 5 females and 2 males)
  are selected.

  A note on singletons: for all per-individual measures, unconnected individuals
  within the pedigree are automatically excluded from the calculations as no
  kinship based statistic can be estimated for them since they do, by definition,
  not share kinship with any other individual in the pedigree.

}
\section{Familial aggregation methods}{

  \describe{

    \item{binomialTest}{
      Evaluate whether the number of affected in a trait are higher than
      expected by chance using a simple binomial test. In contrast to
      most other methods presented here, this does not use the kinship
      between affected individuals, but simply performs a binomial test
      for each family considering the numbers of affected within the
      family, the size of the family and the global probability of being
      affected. The latter is by default calculated on the data set
      (ratio between the total number of affected in the pedigree and
      the total number of phenotyped individuals), can however also be
      specified with the \code{prob} argument.

      The test is performed using the \code{\link{binom.test}}.

      The function returns a \code{\linkS4class{FABinTestResults}} object.
    }
    
    \item{familialIncidenceRate}{
      Calculate the familial incidence rate (FIR, or FR) as defined in [Kerber
      1995], formula (3). The FIR is an estimate for the risk per
      gene-time for each individual for a certain disease (trait) given
      the disease experience in the cohort. The measure considers the
      kinship of each individual with any affected individual in the
      pedigree and the time at risk for each individual.

      Internally, the function first excludes individuals from the test
      which have a missing value (\code{NA}) either in the argument
      \code{trait} or in the argument \code{timeAtRisk}. Next, the thus
      reduced pedigree, is further cleaned by removing all resulting
      singletons (i.e. individuals that do not share kinship with any
      other individual in the above reduced data set).

      The method returns a vector with the FIR value for each
      individual. Individuals that were excluded from the test as
      described above habe an FIR of \code{NA}.
    }

    \item{familialIncidenceRateTest}{
      Calculates the familial incidence rate for each individual and in
      addition assesses the significance of these based on Monte Carlo
      simulations. See \code{\linkS4class{FAIncidenceRateResults}} for
      more details.

      The method returns a \code{\linkS4class{FAIncidenceRateResults}}
      object.
    }

    \item{fsir}{
      Calculate the familial standardized incidence rate (FSIR) as
      defined in [Kerber, 1995], formula (4). The FSIR weights the
      disease status of relatives based on their degree of relatedness
      with the proband [Kerber, 1995]. Formally, the FSIR is defined as
      the standardized incidence ratio (SIR) or standardized morality
      ratio in epidemiology, i.e. as the ratio between the observed and
      expected number of cases, only that both are in addition also
      weighted by the degree of relatedness (i.e. kinship value) between
      individuals in the pedigree.

      Similar to \code{familialIncidenceRate}, the function excludes
      individuals with missing values in any of the arguments
      \code{trait}, \code{timeInStrata} (and optionally \code{strata})
      and all individuals that do not share any kinship with any other
      individual in the pedigree after removing the above individuals.

      The method returns a vector with the FSIR value for each
      individual. Individuals excluded as above describe have a FSIR value
      of \code{NA}.
    }

    \item{fsirTest}{
      Calculates the familial standardized incidence rate (FSIR) for
      each individual and in addition assesses the significance of these
      based on Monte Carlo simulations. See
      \code{\linkS4class{FAStdIncidenceRateResults}} for more details.

      The method returns a \code{\linkS4class{FAStdIncidenceRateResults}}
      object.
    }

    \item{genealogicalIndexTest}{
      Performs the \emph{genealogical index} analysis from [Hill
      1980] (also known as the \emph{genealogical index of
      familiality} or \code{genetic index of familiality}) to identify
      familial clustering of traits (e.g. cancers etc).

      This test calculates the mean kinship among affected individuals in a
      pedigree along with mean kinships of equal sized random control
      sets drawn form the pedigree. The distribution of average kinship
      values among these random sets is used to estimate the probability
      that the observed mean kinship between the affected individuals is due
      to chance.
      The \code{controlSetMethod} argument allows to specify the method
      to define sets of matched control individuals in a pedigree or
      family.

      Note that by default singletons (i.e. unconnected individuals in
      the pedigree) are removed from the pedigree prior the
      analysis. Set \code{rm.singletons=FALSE} if you do not want them
      to be removed.

      The method can also be performed separately for each family within
      the larger pedigree (\code{perFamilyTest=TRUE} to evaluate the
      presence of clustered affected within each family).
      In this case it is also possible to
      use \code{controlSetMethod="getGenerationMatched"} or
      \code{controlSetMethod="getGenerationSexMatched"}, which allows to
      draw random control samples from the same generation(s) than the
      affected are.

      Stratified random sampling can be performed with the \code{strata}
      argument. See details for more information.

      The function returns a \code{\linkS4class{FAGenIndexResults}} object.
    }

    \item{kinshipGroupTest}{
      Performs a familial aggregation test on a subset of a family. The
      idea behind this test is to narrow down the set of controls for
      each affected individual by considering only individuals that are
      as closely related as the most distant affected individual. This
      strategy incorporates more the family structure of the cases and
      is meant to be an alternative to the \code{kinshipSumTest}
      method.

      Initially, for an affected individual \emph{i} a group \emph{C(i)}
      is created that contains all individuals that share kinship as far as
      the most distantly related affected individual. This cluster can be
      interpreted as a circle that is centered at individual \emph{i} with
      radius equal to the most distantly related case. Therefore, the
      cluster defines a narrowed, individual-specific set of individuals
      in which the phenotype is assumed to have been passed on.
      Groups consisting of the same set of affected individuals are
      reduced to a single group (i.e. the group with the smallest total
      number of individuals).

      Based on this definition of groups \emph{C(i)}, we compute two
      statistics by performing Monte Carlo simulations (which optionally
      allow to perform stratified random sampling). During each
      simulation step affected cases are randomly sampled from the
      population.

      1. The ratio test counts per group \emph{C(i)} the number of times
      we observe a higher number of affected individuals in the simulation
      than in the observed case. Dividing this number by the number of
      simulation steps yields immediately the p-value that describes the
      event to observe by chance a higher number of affected individuals
      than in the given case.

      2. The kinship test addresses the degree of relatedness within the
      simulated set by a counting method where we count the number of times
      in a simulation step there is a pair of affected individuals that are
      more closely related than in the observed group \emph{C(i)}. In case
      the closest degree of relatedness is equal in both the simulation
      step and the observed case, we look at the number of pairs found in
      both and count it if this number is higher in the simulation
      step. Again, dividing this count by the number of simulation steps
      readily yields a p-value.

      See also the method \code{runSimulation} for
      \code{\linkS4class{FAKinGroupResults}}.

      The function returns a \code{\linkS4class{FAKinGroupResults}} object.
    }

    \item{kinshipSumTest}{
      Performs a test for familial aggregation based on the sum of kinship
      values between affected cases. This test highlights individuals that
      exhibit a higher than chance relationship to other affected
      individuals, therefore highlighting individuals within families
      aggregating the phenotype. To achieve this, for each affected
      individual the sum of kinship values to all other affected cases is
      computed. In a Monte Carlo simulation this is repeated with the same
      number of cases (and optionally stratified with the \code{strata}
      argument), and the resulting background distribution is used to
      compute p-values for the kinship sums obtained from the observed
      cases. See also the method \code{runSimulation} for
      \code{\linkS4class{FAKinSumResults}}.

      The function returns a \code{\linkS4class{FAKinSumResults}} object.
    }

    \item{probabilityTest}{
      DEPRECATED: this test will be removed in Bioconductor version 3.8
      due to problems and incompatibilities of the \code{gap} package on
      MS Windows systems.
      
      This is only a convience method that calls the \code{gap} package's
      method \code{pfc.sim} to compute probabilities of familial clustering
      of phenotypes [Yu and Zelterman (2002)]. One drawback of that method
      is that it is limited to families with at most 22 individuals. Thus,
      pedigrees need to be split with specialized software such as Jenti
      [Falchi and Fuchsberger ea. (2008)], which within large families
      define cliques that can then be used as input to this algorithm.

      See also method \code{runSimulation} for
      \code{\linkS4class{FAProbResults}}.

      The function returns a \code{\linkS4class{FAProbResults}} object.
    }


  }
}
\section{Utility functions}{

  \describe{

    \item{factor2matrix}{
      Converts a factor into a matrix with columns corresponding to the
      levels and values (cell row i, column j) being either 0 or 1
      depending on whether the ith factor was of the level j. See
      examples below for in or \code{FAStdIncidenceRateResults}.
    }

    \item{estimateTimeAtRisk}{
      Function to calculate the time at risk based on the start date of
      the study or the birth date of an individual (\code{startDate})
      and the study's end date (\code{endDate}), the date of an
      incidence (e.g. date of diagnosis of a cancer
      \code{incidenceDate}) or the death of the individual
      (\code{deathDate}). The time at risk for each individual is
      calculated as the minimal time period between \code{startDate} and
      any of \code{endDate}, \code{incidenceDate} or
      \code{deathDate}. Thus it is also possible to provide just the
      \code{endDate} along with the \code{startDate}, in which case the
      \code{endDate} should be the earliest time point of: end date of
      the study, incidence date or date of death.

      For affected individuals (those for which either an incidence date
      is provided or the value in the optional argument \code{affected}
      is \code{TRUE} or bigger than 0), by default half of the time unit
      is subtracted. For example, a individual that has an incidence
      after 2 days is 1.5 days at risk. The proportion of the time unit
      to subtract can be specified with the argument
      \code{incidenceSubtract}.

      The function returns a numeric vector with the time at risk in
      days.
    }

    \item{sliceAge}{
      Generates a matrix with columns corresponding to age slices/strata
      defined by argument \code{slices} and rows to individuals. Each
      cell in a row represents the time spent by the individual in the
      age slice/strata. See example below.
    }

  }
}

\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}

\author{
  Johannes Rainer, Daniel Taliun, Christian Weichenberger.
}
\references{
  Rainer J, Talliun D, D'Elia Y, Domingues FS and Weichenberger CX
  (2016) FamAgg: an R package to evaluate familial aggregation of traits
  in large pedigrees. \emph{Bioinformatics}.
  
  Hill, J.R. (1980) A survey of cancer sites by kinship in the Utah Mormon
  population. In Cairns J, Lyon JL, Skolnick M (eds):
  \emph{Cancer Incidence in Defined Populations. Banbury Report 4.}
  Cold Spring Harbor, NY: Cold Spring Harbor Laboratory Press, pp 299--318.

  Kerber, R.A. (1995) Method for calculating risk associated with family
  history of a disease. \emph{Genet Epidemiol}, pp 291--301.

  Yu, C. and Zelterman, D. (2002) Statistical inference for familial disease
  clusters. \emph{Biometrics}, pp 481--491

  Falchi, M. and Fuchsberger, C. (2008) Jenti: an efficient tool for mining
  complex inbred genealogies. \emph{Bioinformatics}, pp 724--726
}

\seealso{
  \code{\link[kinship2]{pedigree}},
  \code{\linkS4class{FAData}},
  \code{\linkS4class{FAProbResults}},
  \code{\linkS4class{FAKinGroupResults}},
  \code{\linkS4class{FAKinSumResults}},
  \code{\linkS4class{FAIncidenceRateResults}}
}
\examples{
##########################
##
##  Defining a small pedigree
##
## load the Minnesota Breast Cancer record and subset to the
## first families.
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid==4 | minnbreast$famid==5 |
                    minnbreast$famid==14 | minnbreast$famid==8, ]
mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
## renaming column names
colnames(mbped) <- c("family", "id", "father", "mother", "sex")
## create the FAData object
fad <- FAData(pedigree=mbped)

## We specify the cancer trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

##########################
##
## Familial Incidence Rate
##
## Calculate the FR for each individual given the affected status of
## each individual in trait cancer and the time at risk for each
## participant. We use column "endage" in the minnbreast data.frame
## that specifies the age at the last follow-up or incident cancer as a
##rather impresice estimate for time at risk.
fr <- familialIncidenceRate(fad, trait=tcancer, timeAtRisk=mbsub$endage)

## Plot the distribution of familial rates
plot(density(fr, na.rm=TRUE))

## Perform in addition Monte Carlo simulations to assess the significance
## for the familial incidence rates.
frRes <- familialIncidenceRateTest(fad, trait=tcancer,
                                   timeAtRisk=mbsub$endage,
                                   nsim=500)
head(result(frRes))


##########################
##
## Familial Standardized Incidence Rate:
## Please see examples of FAStdIncidenceRateResults.



##########################
##
##  Perform familial aggregation analyses using the genealogical index
##
gi <- genealogicalIndexTest(fad, trait=tcancer, traitName="cancer",
                            nsim=500)
result(gi)
## A significant clustering of cancer cases was identified in the
## analyzed pedigree.

## Plotting the observed mean kinship and the distribution of mean kinship
## from the random sampling.
plotRes(gi)


##########################
##
##  Perform familial aggregation analysis using the kinship sum test
##
kcr <- kinshipSumTest(fad, trait=tcancer, traitName="cancer",
                      nsim=500)
kcr
head(result(kcr))


##########################
##
##  Perform familial aggregation analysis using the kinship group test,
##  stratifying by sex
##
kr <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer",
                       nsim=500, strata=fad$sex)
kr
head(result(kr))



##########################
##
##  Estimate the time at risk given
##
## Define some birth dates and incidence dates and end date of study
bdates <- c("2012-04-17", "2014-05-29", "1999-12-31", "2002-10-10")
idates <- c(NA, NA, "2007-07-13", "2013-12-23")
edates <- rep("2015-09-15", 4)

## Estimate the time at risk. The time period is returned in days.
riskDays <- estimateTimeAtRisk(startDate=bdates, incidenceDate=idates, endDate=edates)
riskDays

##########################
##
##  Define the time spent in an age stratum given the indivduals'
##  age at incidence or end of study.
head(mbsub$endage)
## We "slice" the age in specified intervals/slices
stratAge <- sliceAge(mbsub$endage, slices=c(0, 40, 60, Inf))
head(stratAge)

## The first column lists the number of years spent in the first age
## stratum (0 < age <= 40) and the second in the second stratum
## (40 < age <= Inf)

## We could also stratify the disk days from above in per year strata.
sliceAge(riskDays/365, slices=c(0, 2.5, 5, 10, 20))


##########################
##
##  Simple example for factor2matrix: generate a matrix for factor $sex
head(factor2matrix(fad$sex))




}
\keyword{classes}
