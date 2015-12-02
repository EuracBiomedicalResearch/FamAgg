\name{FAKinGroupResults-class}
\Rdversion{1.1}
\docType{class}
\alias{FAKinGroupResults-class}
\alias{FAKinGroupResults}
\alias{[,FAKinGroupResults,ANY-method}
\alias{[,FAKinGroupResults,ANY,ANY,ANY-method}
\alias{affectedKinshipGroups}
\alias{affectedKinshipGroups,FAKinGroupResults-method}
\alias{buildPed,FAKinGroupResults-method}
%\alias{independentGroupCount}
%\alias{independentGroupCount,FAKinGroupResults-method}
\alias{plotPed,FAKinGroupResults-method}
\alias{plotRes,FAKinGroupResults-method}
\alias{result,FAKinGroupResults-method}
\alias{runSimulation,FAKinGroupResults-method}
\alias{shareKinship,FAKinGroupResults-method}
\alias{trait<-,FAKinGroupResults-method}


\title{Kinship group test}

\description{
  The \code{FAKinGroupResults} object contains the results from the
  \emph{kinship test}. This test performs a familial aggregation
  analysis on a subset of individuals within a family. Two actual tests
  are conducted, a \emph{ratio test} that evaluates whether the number
  of affected individuals within the group is higher than expected by
  chance, and a \emph{kinship test} that compares the largest kinship
  value between affected in the group to the one between randomly
  sampled individuals.

  For more details see \code{\link{kinshipGroupTest}}.
}
\section{Objects from the Class}{
  \code{FAKinGroupResults} objects are created by the
  \code{\link{kinshipGroupTest}} method for \code{\linkS4class{FAData}}
  objects.
}
\usage{

\S4method{affectedKinshipGroups}{FAKinGroupResults}(object)

\S4method{buildPed}{FAKinGroupResults}(object, id=NULL, max.generations.up=3,
                                     max.generations.down=16, prune=FALSE)

% \S4method{independentGroupCount}{FAKinGroupResults}(object, method="jo")

\S4method{plotPed}{FAKinGroupResults}(object, id=NULL, family=NULL,
                                    filename=NULL, device="plot", ...)

\S4method{plotRes}{FAKinGroupResults}(object, id=NULL, family=NULL,
                                    addLegend=TRUE, type="density", ...)

\S4method{result}{FAKinGroupResults}(object, method="BH")

\S4method{runSimulation}{FAKinGroupResults}(object, nsim=50000, strata=NULL)

\S4method{shareKinship}{FAKinGroupResults}(object, id=NULL)

\S4method{trait}{FAKinGroupResults}(object) <- value

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
    For \code{buildPed}, \code{plotPed}, \code{plotRes} and
    \code{shareKinship}: the id of the kinship group (i.e. one of the
    ids in column \code{"group_id"} of the result table \code{result(object)}).

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
    %% The method to calculate the number of independent groups of affected
    %% individuals in the pedigree. Method \code{"jo"} returns the number
    %% of independent groups, i.e. groups that do not share a single
    %% affected individual. Method \code{"daniel"} uses the original code
    %% from Daniel Taliun.
  }

  \item{nsim}{
    Number of simulations.
  }

  \item{object}{
    The \code{FAKinGroupResults} object.
  }

  \item{prune}{
    For \code{buildPed}: whether the full pedigree should be returned
    (\code{prune=FALSE}) or the pedigree should be reduced to a smaller
    pedigree containing only individuals in the kinship group
    (\code{prune=TRUE}); see details for more information. Note: the
    \code{plotPed} method does also support this parameter.
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
      The list containing the result of the simulation. Element
      \code{"pvalueKinship"} contains the p-values from the kinship test,
      \code{"pvalueRatio"} the p-values from the ratio test,
      \code{"expDensity"} and \code{"expHist"} \code{density} and
      \code{hist} objects representing the background distribution from
      the Monte Carlo simulation.
      This slot should not be directly accessed, use the \code{result}
      method to extract result information.
    }
    \item{affectedKinshipGroups}{
      A list of lists, each element representing one kinship group,
      \code{$aff} the ids of all affected persons in that group sharing
      kinship with the indivudual (being the id of the group),
      \code{$phe} ids of all phenotyped individuals in the pedigree with
      kinship to the individual up to a kinship being smaller or equal
      to the smallest kinship of the indivudal with any other affected
      in that group. \code{$kinfreq}: a table with the frequency
      (counts) of kinship values (smaller 0.5). \code{$meankin}: the
      mean kinship in that group.
    }
  }
}
\section{Methods and Functions}{
  \describe{

    \item{affectedKinshipGroups}{
      Get groups of affected individuals in the pedigree along with
      all phenotyped individuals with kinship larger or equal to the
      smallest kinship between affected individuals in the
      group. Returns a \code{list} of lists with elements \code{aff},
      \code{phe}, \code{kinfreq} and \code{meankin} for each list item:
      \code{aff}: character vector with the ids of all affected in the
      group, \code{phe}: character vector with the ids of all phenotyped
      in the group, \code{kinfreq}: a \code{table} with the frequency
      (counts) of kinship values (self-self kinships removed; the names of the
      table correnspond to the kinship values ordered increasingly), and
      \code{meankin}: the mean kinship value in the group. The names of
      the list correspond to the id of the affected individual for which
      the affected kinship group was determined.
    }

    \item{buildPed}{
      Builds the pedigree for the submitted id (which represents the id
      of the group, i.e. one of
      \code{names(affectedKinshipGroups(object))}, respectively ids in
      column \code{"group_id"} of \code{result(object)}). Building a
      pedigree by submitting the family id is not supported. See details
      below for more information. The resulting pedigree is returned as
      a \code{data.frame}.
      %% Builds the pedigree for the submitted id (which represents the id
      %% of the group, i.e. one of
      %% \code{names(affectedKinshipGroups(object))}). In contrast to the
      %% same method for \code{\linkS4class{FAData}} objects, this method
      %% includes only individuals in the pedigree that are in the selected
      %% group of individuals with eventually included missing mates for
      %% family founders.
    }

    %% \item{independentGroupCount}{
    %%   Returns the number of independent groups of affected individuals.
    %% }

    \item{plotPed}{
      Plots a pedigree for one of the affected kinship groups in the
      simulation results. The ids of all individuals of the affected
      kinship group are highlighted in the plot in red.
      See \code{\link{plotPed}} for more details on
      the plotting and details below for additional settings.
    }

    \item{plotRes}{
      Plots the distribution of counts of randomly sampled affected
      counts within the kinship group along with the
      actually observed kinship sum.
    }

    \item{result}{
      Returns the result from the simulation as a \code{data.frame} with
      columns:
      \code{"trait_name"}: the name of the trait.
      \code{"total_phenotyped"}: the total number of phenotyped
      individuals in the trait.
      \code{"total_affected"}: the total number of affected individuals
      in the trait.
      \code{"phenotyped"}: the number of phenotyped individuals in all
      analyzed kinship groups.
      \code{"affected"}: the number of affected individuals in all
      analyzed kinship groups.
      \code{"group_id"}: the id for the kinship group (represents the id
      of one of the affected individuals in the group).
      \code{"family"}: the family id of the affected/group.
      \code{"group_phenotyped"}: the number of phenotyped individuals in the
      current kinship group.
      \code{"group_affected"}: the number of affected individuals in the
      current kinship group.
      \code{"ratio_pvalue"}: the p-value from the ratio test.
      \code{"ratio_padj"}: the p-value from the ratio test adjusted for
      multiple hypothesis testing using the method specified with
      argument \code{method}.
      \code{"mean_kinship"}: the mean kinship value between all
      individuals in the current kinship group.
      \code{"kinship_pvalue"}: the p-value from the kinship test.
      \code{"kinship_padj"}: the p-value from the kinship test adjusted
      for multiple hypothesis testing using the method specified with
      argument \code{method}.

      The \code{data.frame} is sorted by column \code{"ratio_pvalue"},
      its rownames correspond to the \code{"group_id"}.
    }

    \item{runSimulation}{
      Performs the simulation analysis based on the pedigree and trait
      information stored in the object. Optionally allows to perform
      stratified sampling. Returns a \code{FAKinGroupResults}
      with the results of the simulation.
    }

    \item{shareKinship}{
      Returns a character vector with ids of all individuals that share
      kinship with any of the individuals in the kinship group
      identified by the argument \code{id}.
    }

    \item{trait<-}{
      Set the trait information. This method will reset all simulation
      results saved in the \code{sim} slot.
    }

  }
}
\details{
  Calling the \code{runSimulation} method on a \code{FAKinGroupResults}
  object is the same as calling the \code{\link{kinshipGroupTest}} on a
  \code{\linkS4class{FAData}} object. In the first case the
  simulation is performed using the trait information
  data stored internally in the object, while in the latter case the
  trait information have to be submitted to the function
  call.

  A call to the setter methods \code{trait<-} resets any simulation
  results present in the \code{sim} slot, thus, the object can be
  re-used to perform a simulation analysis using the new trait data.

  The \code{buildPed} method returns by default the full pedigree (all
  ancestors and all children) up to the maximal number of
  generations. By setting \code{prune=TRUE} the method restricts the
  pedigree to all individuals with a kinship >= the minimal kinship
  between the individual (with the id equal to the group id) and any
  other affected individual in its pedigree.

  The \code{plotPed} method allows to plot the pedigree for a kinship
  group. This pedigree consists of the full pedigree (all ancestors and
  children) of all individuals in the kinship group. Similar to the
  \code{buildPed} method the pedigree can be restricted to the kinship
  group (and eventual missing parents etc) by setting \code{prune=TRUE}.
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
\seealso{
  \code{\linkS4class{FAData}},
  \code{\link{kinship}},
  \code{\link{trait}},
  \code{\link{probabilityTest}},
  \code{\link{kinshipSumTest}},
  \code{\link{genealogicalIndexTest}},
  \code{\link{familialIncidenceRateTest}},
  \code{\link{fsirTest}},
  \code{\link{buildPed}},
  \code{\link{plotPed}},
  \code{\link{switchPlotfun}}
}
\examples{
##########################
##
##  Perform the simulation analysis
##
## Load the test data.
data(minnbreast)

## Subset to some families and generate the pedigree data.frame
mbsub <- minnbreast[minnbreast$famid == 165 | minnbreast$famid == 432, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- c("family", "id", "father", "mother", "sex")

## Generate the FAData.
fad <- FAData(pedigree=PedDf)

## Specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## Perform the kinship group test.
far <- kinshipGroupTest(fad, trait=tcancer, traitName="cancer", nsim=1000)
res <- result(far)
head(res)

## Plot the pedigree for the most significant kinship group
plotPed(far, id=res[1, "group_id"])

## The full pedigree for this affected individual and its kinship group is
## large:
nrow(buildPed(far, id=res[1, "group_id"]))

## We can however restrict it to a reduced pedigree containing only the
## kinship group and all individuals with a kinship >= the smallest kinship
## between the individual and any other affected individual in the pedigree:
nrow(buildPed(far, id=res[1, "group_id"], prune=TRUE))

## By specifying prune=TRUE we can restrict the pedigree plot to these
## individuals
plotPed(far, id=res[1, "group_id"], prune=TRUE)

## Get the ids of all individuals sharing kinship with any of the inddividuals
## in that kinship group:
shareKinship(far, id=res[1, "group_id"])

## Plot the simulation analysis result for the ratio test.
plotRes(far, id=res[1, "group_id"], type="hist")

}
\keyword{classes}

