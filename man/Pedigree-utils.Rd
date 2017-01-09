\name{PedigreeUtils}
\Rdversion{1.1}
\alias{PedigreeUtils}
\alias{connectedSubgraph}
\alias{countGenerations}
\alias{countGenerations,FAData-method}
\alias{countGenerations,data.frame-method}
\alias{countGenerations,pedigreeList-method}
\alias{countGenerations,pedigree-method}
%% estimateGenerations
\alias{estimateGenerations}
\alias{estimateGenerations,FAData-method}
\alias{estimateGenerations,data.frame-method}
\alias{estimateGenerations,pedigreeList-method}
\alias{estimateGenerations,pedigree-method}
%% findFounders
\alias{findFounders}
\alias{findFounders,FAData-method}
\alias{findFounders,data.frame-method}
\alias{findFounders,pedigreeList-method}
\alias{findFounders,pedigree-method}
%% generationsFrom
\alias{generationsFrom}
\alias{generationsFrom,FAData-method}
\alias{generationsFrom,data.frame-method}
\alias{generationsFrom,pedigreeList-method}
\alias{generationsFrom,pedigree-method}
%% getAncestors
\alias{getAncestors}
\alias{getAncestors,FAData-method}
\alias{getAncestors,data.frame-method}
\alias{getAncestors,pedigreeList-method}
\alias{getAncestors,pedigree-method}
%% getChildren
\alias{getChildren}
\alias{getChildren,FAData-method}
\alias{getChildren,data.frame-method}
\alias{getChildren,pedigreeList-method}
\alias{getChildren,pedigree-method}
%% getCommonAncestor
\alias{getCommonAncestor}
\alias{getCommonAncestor,FAData-method}
\alias{getCommonAncestor,data.frame-method}
\alias{getCommonAncestor,pedigreeList-method}
\alias{getCommonAncestor,pedigree-method}
%% getFounders
\alias{getFounders}
\alias{getFounders,FAData-method}
\alias{getFounders,data.frame-method}
%% getMissingMate
\alias{getMissingMate}
\alias{getMissingMate,FAData-method}
\alias{getMissingMate,data.frame-method}
\alias{getMissingMate,pedigreeList-method}
\alias{getMissingMate,pedigree-method}
%% getSiblings
\alias{getSiblings}
\alias{getSiblings,FAData-method}
\alias{getSiblings,data.frame-method}
\alias{getSiblings,pedigreeList-method}
\alias{getSiblings,pedigree-method}
%% getSingletons
\alias{getSingletons}
\alias{getSingletons,FAData-method}
\alias{getSingletons,data.frame-method}
%% removeSingletons
\alias{removeSingletons}
\alias{removeSingletons,FAData-method}
\alias{removeSingletons,data.frame-method}

\alias{ped2graph}
\alias{subPedigree}
\alias{cliques,missing-method}

\alias{shareKinship}
\alias{shareKinship,FAData-method}


\title{Basic pedigree utilities}

\description{
  Utility functions to access, modify or subset pedigrees. Most of these
  functions can be applied to simple \code{data.frame} in pedigree
  format or \code{pedigree} or \code{pedigreeList} objects defined in
  the \code{kinship2} package.
}
\usage{

\S4method{cliques}{missing}(object, ...)

connectedSubgraph(graph, nodes, mode="all", all.nodes=TRUE, ifnotfound)

\S4method{countGenerations}{FAData}(object, id=NULL, direction="down", ...)

\S4method{estimateGenerations}{FAData}(object, family=NULL, ...)

\S4method{findFounders}{FAData}(object, family=NULL, id = NULL, ...)

\S4method{generationsFrom}{FAData}(object, id=NULL, ...)

\S4method{getAncestors}{FAData}(object, id=NULL, max.generations=3, ...)

\S4method{getChildren}{FAData}(object, id=NULL, max.generations=16, ...)

\S4method{getCommonAncestor}{FAData}(object, id, method="min.dist")

\S4method{getFounders}{FAData}(object, ...)

\S4method{getMissingMate}{FAData}(object, id=NULL, ...)

\S4method{getSiblings}{FAData}(object, id=NULL, ...)

\S4method{getSingletons}{FAData}(object, ...)

ped2graph(ped)

\S4method{removeSingletons}{FAData}(object, ...)

\S4method{removeSingletons}{data.frame}(object, ...)

subPedigree(ped, id=NULL, all=TRUE)

\S4method{shareKinship}{FAData}(object, id)

}
\arguments{

  (in alphabetic order)

  \item{all}{
    For \code{subPedigree}: if all individuals have to be present in the
    sub-pedigree.
  }

  \item{all.nodes}{
    For \code{connectedSubgraph}: if all nodes have to be present in the
    resulting graph, or only those that are connected with each other.
  }

  \item{direction}{
    For \code{countGenerations}: whether the number of ancestor (\code{"up"})
    generations or offspring (\code{"down"}) generation should be counted.
  }

  \item{family}{
    A character or numeric representing the family id.
    For \code{doFindFounders}: the id of the family in the pedigree for
    which the founders should be identified. Uses the first family in
    the pedigree if not specified.
    For \code{estimateGenerations}: optional id of the family if
    generation numbers should only be calculated for one
    family. Otherwiseif, the generations are estimated for all families
    (separately) in the object.
  }

  \item{graph}{
    An \code{igraph} graph object.
  }

  \item{id}{
    A character or numeric vector length 1 or longer specifying the
    id(s) of the individual(s).
    For \code{generationsFrom} and \code{findFounders} only a single id
    should be submitted.
  }

  \item{ifnotfound}{
    For \code{connectedSubgraph}: if not defined, the function throws an
    error if no subgraph can be specified. If defined, its value is
    returned if no subgraph was found.
  }

  \item{max.generations}{
    For \code{getAncestors} and \code{getChildren}: the maximal
    number of ancestor or offspring generations that should be returned.
  }

  \item{method}{
    For \code{getCommonAncestor}: the method by which the closest
    common ancestor sould be identified. Either \code{"min.dist"}
    (ancestor with the minimal distance to any of the individuals) or
    \code{"smallest.mean.dist"} (ancestor with the smallest mean
    distance to any of the individuals).
  }

  \item{mode}{
    For \code{connectedSubgraph}: either \code{"all"}, \code{"in"},
    \code{"out"} specifying how distances and paths between individual
    nodes should be determined. See help for function
    \code{shortest_paths} in package \code{igraph} for more details.
  }

  \item{nodes}{
    For \code{connectedSubgraph}: A character vector of node (vertex) names
    for which the subgraph should be defined.
  }

  \item{object}{
    For \code{cliques}: passed to the \code{cliques} function
    from the \code{igraph} package.
    For all other methods: either a \code{\linkS4class{FAData}} object
    (or any object inheriting from it), a \code{data.frame},
    \code{pedigree} or \code{pedigreeList} objects (the latter being
    defined in the \code{kinship2} package).
  }

  \item{ped}{
    Either a \code{data.frame} or a \code{pedigree} object specifying
    the pedigree. If a \code{data.frame} is submitted, the columns
    \code{id}, \code{family}, \code{father}, \code{mother} and
    \code{sex} are required.
  }

  \item{...}{
    For \code{cliques}: additional arguments passed to the
    \code{cliques} function from the \code{igraph} package.
  }

}
\section{Basic pedigree utilities}{

  \describe{

    \item{countGenerations}{
      Count the generations up- or down the pedigree for the specified
      individual(s), i.e. determine the number of ancestor or offspring
      generations defined in the pedigree for the specified
      individual(s).
      Returns a named numeric vector, names corresponding to the
      individual's id, with the number of generations for each specified
      individual.
    }

    \item{findFounders}{
      Identifies the founder couple with the largest number of offspring
      generations in the pedigree. The provided pedigree
      object/data.frame can contain pedigrees of multiple families,
      thus, to identify the founder pair for a family its ID can be
      provided with the \code{family} parameter. Alternatively, the ID
      of an individual can be specified, in which case the founder pair
      of the (full) pedigree of the specified individual is identified.
      If two or more couples have the same, largest number of
      offspring generations, the first couple is selected.
      Returns a character vector of length 2 with the ids of the founder
      individuals.
    }

    \item{getFounders}{
      Returns the ids of all founders in the pedigree. A founder is an
      individual from which neither father nor mother is known in the
      pedigree.
    }

    \item{getSingletons}{
      Returns the ids of all singletons, i.e. individuals in the
      pedigree that are not connected to any other individual (have no
      parents in the pedigree and no children).
    }

    \item{getAncestors}{
      Identify and return the ids of ancestor generations (up to
      \code{max.generations}) for the specified individual(s).
    }

    \item{getChildren}{
      Identify and return the ids of offspring generations (up to
      \code{max.generations}) for the specified individual(s).
    }

    \item{getCommonAncestor}{
      Finds the closest common ancestor between specified individuals (2
      or more ids are required).
      Returns a character vector with the ids of the ancestors or \code{NA} if
      no common ancestor was found.
    }

    \item{getMissingMate}{
      The function evaluates if in the sub-pedigree defined by the
      specified ids one or more mates (spouse) are
      missing and if so it returns their ids.
    }

    \item{getSiblings}{
      Get siblings for the specified id(s). Returns their ids as
      character, or numeric vector.
    }

    \item{removeSingletons}{
      Removes all unconnected individuals (i.e. singletons) from the
      pedigree. Returns a \code{data.frame} with the pedigree cleaned
      from all singletons. Note that, due to internal sanitizing,
      columns \code{"father"} and \code{"mother"} in the resulting
      \code{data.frame} have a \code{NA} for individuals for which the
      father or mother is not known in the pedigree.
    }

    \item{subPedigree}{
      Finds the smallest pedigree containing all specified individuals.
      Depending on the input, a \code{data.frame}, \code{pedigree}
      or \code{pedigreeList}.
    }

  }
}

\section{Advanced pedigree methods}{
  \describe{

    \item{estimateGenerations}{
      Estimates generation levels/numbers for each, or only one, family
      in the object. Generation numbers are always relative to the
      founder couple (defined by \code{\link{findFounders}}).
      Returns (always) a named list of generation numbers. The names of
      the list represent the family id, the names of the numeric vector
      of generations the id of the individuals in the family.
    }

    \item{generationsFrom}{
      Determine generations starting from the specified
      individual. Siblings including their mates and all other in the
      same generation () are assigned generation 0, ancestor generations
      (all their parents etc) negative generation numbers, decreasing with
      ancestor level and their offspring positive numbers, increasing
      with each generation.
      Generations are only estimated within the family of the
      individual, also, if the pedigree consists of un-connected
      sub-pedigree, generation numbers will only be calculated for the
      sub-pedigree containing the specified individual.
      The function returns a named numeric vector of generation numbers,
      the names corresponding to the ids of the individuals in the
      specified individual's family. Not connected individuals in the
      family get a \code{NA} generation number.
    }

    \item{shareKinship}{
      Finds all related individuals (individuals sharing kinship with
      the individual) for the specified individual(s) in the pedigree
      and returns their ids as a character vector.
    }
  }
}

\section{Graph theory related functions}{
  \describe{
    \item{cliques}{
      Wrapper method passing all arguments to the \code{cliques}
      function from the \code{igraph} package.
    }

    \item{connectedSubgraph}{
      Finds the (eventually smallest) connected subgraph of all
      specified nodes.
      Returns an \code{igraph} object representing the subgraph of the
      specified nodes.
    }

    \item{ped2graph}{
      Transforms the pedigree into a (directed) graph with the direction
      of the edges being always from parent to child.
      An \code{igraph} object.
    }
  }
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer.
}
\seealso{
  \code{\link[kinship2]{pedigree}},
  \code{\linkS4class{FAData}},
  \code{\linkS4class{FAProbResults}},
  \code{\linkS4class{FAKinGroupResults}},
  \code{\linkS4class{FAKinSumResults}},
  \code{\link{PedigreeAnalysis}}
}
\examples{
##########################
##
##  Defining a small pedigree
##
## load the Minnesota Breast Cancer record and subset to the
## first families.
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid==4 | minnbreast$famid==5, ]
mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
## renaming column names
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

## plot the pedigree for family 4 to get an overview.
switchPlotfun(method="ks2paint")
fam4 <- mbped[mbped$family==4, ]
doPlotPed(individual=fam4$id, father=fam4$father, mother=fam4$mother,
          gender=fam4$sex, device="plot")

## find the closest common ancestor between individuals 23, 3 and 8
getCommonAncestor(fam4, id=c(23, 3, 8))

## create the smallest sub-pedigree for individuals 21, 22 and 25
subPedigree(fam4, id=c(21, 22, 25))
## plot that
fam4sub <- subPedigree(fam4, id=c(21, 22, 25))
doPlotPed(individual=fam4sub$id, father=fam4sub$father, mother=fam4sub$mother,
          gender=fam4sub$sex, device="plot")

#########################
##
##  Basic pedigree utils
##
## Note: the same methods can be applied to a data.frame representing
## a pedigree, or a FAData, pedigree or pedigreeList object.

## Find the founder couple for family 4
findFounders(fam4, family=4)

## Alternatively, find the founders for the pedigree in which ibdividual 20 is a
## member
findFounders(fam4, id = 20)

## Return all founders in the pedigree.
getFounders(fam4)

## Get all founders without children (i.e. singletons).
getSingletons(fam4)

## Clean the pedigree from all singletons
fam4noS <- removeSingletons(fam4)
nrow(fam4)
nrow(fam4noS)

## Count the offspring generations for individual "4"
countGenerations(fam4, id="4")

## Get the ids of all ancestors for that individual
getAncestors(fam4, id="4")

## Get the ids of the children of this individual
getChildren(fam4, id="4", max.generations=1)

## Get the ids of the complete offspring for this individuals
getChildren(fam4, id="4")

## Create a FAData object from the pedigree data.frame
fad <- FAData(fam4)
## get the list of all ids sharing kinship with individuals
## 5 and 9
shareKinship(fad, id=c("5", "9"))

## Count the numbers of generations of ancestors for individual 12
countGenerations(fad, id="12", direction="up")

## Count the numbers of offspring generations for individuals 2 and 29
countGenerations(fad, id=c("2", "29"))

## Get all brothers/sisters for individual 9
getSiblings(fad, id="9")

## Determine generation levels starting from individual "9"
generationsFrom(fad, id="9")

## Estimate generations relative to the founder couple for each
## family in the submitted object, a data.frame in the example below
estimateGenerations(mbped)


#########################
##
##  Graph utilities
##
## Convert the pedigree into a graph
pgraph <- ped2graph(fam4)
plot(pgraph)

## Make a subgraph containing nodes 10, 22, 12 and 14
sgraph <- connectedSubgraph(pgraph, c("10", "22", "12", "14"))
plot(sgraph)

}
\keyword{classes}

