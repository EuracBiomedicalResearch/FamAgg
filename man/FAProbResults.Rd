\name{FAProbResults-class}
\Rdversion{1.1}
\docType{class}
\alias{[,FAProbResults,ANY-method}
\alias{[,FAProbResults,ANY,ANY,ANY-method}
\alias{FAProbResults-class}
\alias{FAProbResults}
\alias{buildPed,FAProbResults-method}
\alias{cliqueAndTrait}
\alias{cliqueAndTrait,FAProbResults-method}
\alias{cliques}
\alias{cliques,FAProbResults-method}
\alias{cliques<-}
\alias{cliques<-,FAProbResults-method}
\alias{plotPed,FAProbResults-method}
\alias{result}
\alias{result,FAProbResults-method}
\alias{runSimulation}
\alias{runSimulation,FAProbResults-method}
\alias{shareKinship,FAProbResults-method}
\alias{trait<-,FAProbResults-method}
\alias{traitByClique}
\alias{traitByClique,FAProbResults-method}


\title{Probability test}

\description{
  The \code{FAProbResults} object contains the results from the
  \emph{probability test}.
  The probability test is only a convience method that calls the \code{gap} package's
  method \code{pfc.sim} to compute probabilities of familial clustering
  of phenotypes [Yu and Zelterman (2002)]. One drawback of that method
  is that it is limited to families with at most 22 individuals. Thus,
  pedigrees need to be split with specialized software such as Jenti
  [Falchi and Fuchsberger ea. (2008)], which within large families
  define cliques that can then be used as input to this algorithm.
}
\section{Objects from the Class}{
  \code{FAProbResults} objects are created by the
  \code{\link{probabilityTest}} method of \code{\linkS4class{FAData}}
  objects.
}
\usage{

\S4method{buildPed}{FAProbResults}(object, id=NULL, max.generations.up=3,
                                  max.generations.down=16, prune=FALSE)

\S4method{cliqueAndTrait}{FAProbResults}(object, na.rm=FALSE)

\S4method{cliques}{FAProbResults}(object, na.rm=FALSE)

\S4method{cliques}{FAProbResults}(object) <- value

\S4method{plotPed}{FAProbResults}(object, id=NULL, family=NULL,
                                 filename=NULL, device="plot", ...)

\S4method{result}{FAProbResults}(object, method="BH")

\S4method{runSimulation}{FAProbResults}(object, nsim=50000)

\S4method{shareKinship}{FAProbResults}(object, id=NULL)

\S4method{trait}{FAProbResults}(object) <- value

\S4method{traitByClique}{FAProbResults}(object)

}
\arguments{

  (in alphabetic order)

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
    For \code{buildPed}, \code{plotPed}, \code{shareKinship}: the id
    (character or numerif) of the clique (i.e. one of the ids in column
    \code{"group_id"} of the result table \code{result(object)}).

    Note: \code{id} can be a \code{numeric} or a
    \code{character}. Numeric ids will be internally converted to
    \code{character}.
  }

  \item{max.generations.down}{
    For \code{buildPed}: the maximal number of generations to look for
    children.
  }

  \item{max.generations.up}{
    For \code{buildPed}: the maximal number of generations to look for
    ancestors.
  }

  \item{method}{
    The multiple hypothesis testing method. All methods supported by
    \code{\link{p.adjust}} are allowed.
  }

  \item{na.rm}{
    Whether \code{NA} elements should be returned or not.
  }

  \item{nsim}{
    Number of simulations.
  }

  \item{object}{
    The \code{FAProbResults} object.
  }

  \item{prune}{
    For \code{buildPed}: whether the full pedigree should be returned
    (\code{prune=FALSE}) or the pedigree should be reduced to the
    individuals in the corresponding clique (\code{prune=TRUE}). Note: the
    \code{plotPed} method does also support this parameter.
  }

  \item{value}{
    For \code{cliques<-}: can be a named numeric, character or factor
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
      The result of the simulation. This slot should not be directly
      accessed, use the \code{result} method to extract result information.
    }

    \item{.cliques}{
      A factor with the assignment of individuals to cliques. This slot
      should not be accessed directly, but the \code{cliques} and
      \code{cliques<-} methods should be used instead that ensure that
      the data is matched to the information in the pedigree.
    }

  }
}
\section{Methods and Functions}{
  \describe{

    \item{cliqueAndTrait}{
      Get a \code{data.frame} with the clique ID and the value from the
      trait for each individual. If \code{na.rm=TRUE} all rows
      (individuals) with either a missing clique ID or trait value are
      removed.
    }

    \item{buildPed}{
      Builds the pedigree for the submitted id (which represents the id
      of the group, i.e. the clique (e.g. column \code{"group_id"} of
      \code{result(object)}). By default the method builds the full
      pedigree for all individuals in the clique but the argument
      \code{prune} allows to reduce it to the individuals of the clique.
      Building a pedigree by submitting the family id is not supported.
      The pedigree is returned as a \code{data.frame}.
    }

    \item{cliques}{
      Returns a \code{factor} vector representing the clique/group
      assignment of the individuals. By default, the ordering
      is the same as pedigree, setting argument \code{na.rm=TRUE}
      removes all \code{NA} values, thus the ordering and length might
      differ (e.g. if some individuals are not part of any clique).
    }

    \item{cliques<-}{
      Setter for the cliques vector. Can be a named numeric vector, character
      vector or factor with the names matching the ids of the
      individuals in the pedigree. The method internally matches and
      re-orders the cliques vector to match the ordering of the ids in
      the pedigree.
      In addition, the function resets eventually stored simulation
      results in the \code{sim} slot.
    }

    \code{plotPed}{
      Plots a pedigree for one of the cliques in the
      simulation results. The ids of all individuals of the selected
      clique are highlighted in red.
      See \code{\link{plotPed}} for more details on
      the plotting and details below for additional settings.
    }

    \item{result}{
      Returns the result from the simulation as a \code{data.frame} with
      columns:
      \code{trait_name}: the name of the trait.
      \code{total_phenotyped}: the total number of phenotyped
      individuals in the trait.
      \code{total_affected}: the total number of affected individuals in
      the trait.
      \code{phenotyped}: the number of phenotyped individuals in
      all group (i.e. cliques specified by the \code{clique} parameter).
      \code{affected}: the number of affected individuals in
      all group (i.e. cliques specified by the \code{clique} parameter).
      \code{group_id}: the id of the group (clique).
      \code{group_phenotyped}: the number of phenotyped individuals in the
      current group.
      \code{group_affected}: the number of affected individuals in the current
      group.
      \code{pvalue}: the p-value from the Monte Carlo simulation.
      \code{padj}: the p-value adjusted for multiple hypothesis testing
      using the method specified with argument \code{method}.
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree, trait and
      clique information stored in the object. Returns a
      \code{FAProbResults} object with the results of the simulation.
    }

    \item{shareKinship}{
      Returns a character vector with the ids of all individuals that
      share kinship with any of the individuals in the clique identified
      by the argument \code{id}.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

    \item{traitByClique}{
      Summarize the trait information by clique. Returns a \code{matrix}
      with the size of the clique (i.e. individuals in the clique with
      available trait information) and count of individuals in the
      clique with a trait value other than zero (column
      \code{affected_size}). The clique ids are used as rownames of the
      \code{matrix}.
    }
  }
}
\details{
  Calling the \code{runSimulation} method on a \code{FAProbResults}
  object is the same as calling the \code{probabilityTest} on a
  \code{\linkS4class{FAData}} object. In the first case the
  simulation is performed using the clique and trait information
  data stored internally in the object, while in the latter case the
  clique and trait information have to be submitted to the function
  call.

  A call to the setter methods \code{trait<-} or \code{cliques<-} resets
  any simulation results present in the \code{sim} slot, thus, the
  object can be re-used to perform a simulation analysis using the new
  trait or clique data.
}
\note{
  Subsetting (using the \code{[} operator) is not supported.
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer, Daniel Taliun
}
\references{
  Yu C & Zelterman D (2002) Statistical inference for familial disease
  clusters. \emph{Biometrics}, pp 481-491

  Falchi M & Fuchsberger C (2008) Jenti: an efficient tool for mining
  complex inbred genealogies. \emph{Bioinformatics}, pp 724-726
}
\seealso{
  \code{\linkS4class{FAData}},
  \code{\link{buildPed}},
  \code{\link{plotPed}},
  \code{\link{trait}},
  \code{\link{probabilityTest}},
  \code{\link{kinshipGroupTest}},
  \code{\link{kinshipSumTest}},
  \code{\link{genealogicalIndexTest}},
  \code{\link{familialIncidenceRateTest}},
  \code{\link{fsirTest}}
}
\examples{

## Build the data based on the internal test data set.
pedFile <- system.file("txt/Large.ped", package="FamAgg")
cliqFile <- system.file("txt/Large-cliques.txt", package="FamAgg")
ped <- read.table(pedFile, sep=";")
cliqDf <- read.table(cliqFile, sep="\t")
cliqs <- cliqDf[, 1]
names(cliqs) <- cliqDf[, 2]

## generate the FAData
fad <- FAData(ped[, 1:5])

## define the trait
affs <- ped[, 6]
affs[affs==2] <- 1
names(affs) <- ped[, 2]

## Performing the test
far <- probabilityTest(fad, trait=affs,
                       cliques=cliqs, nsim=1000)

## display the results:
res <- result(far)
head(res)

## plot the pedigree for the first group
## we are using the kinshiup2 plotting function, since the pedigree is
## quite busy
plotPed(far, id=res[1, "group_id"])

## the full pedigree for this affected individual is:
nrow(buildPed(far, id=res[1, "group_id"]))

## restrict to the individuals in the tested clique:
nrow(buildPed(far, id=res[1, "group_id"], prune=TRUE))
plotPed(far, id=res[1, "group_id"], prune=TRUE)


}
\keyword{classes}

