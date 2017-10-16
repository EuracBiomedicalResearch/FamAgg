\name{plotPed}
%\Rdversion{1.1}
%\docType{class}
\alias{doPlotPed}
\alias{switchPlotfun}

\title{Pedigree plotting functions}

\description{
  Plot a pedigree for a family or an individual.
}
\usage{

doPlotPed(family=NULL, individual=NULL, father=NULL, mother=NULL, gender=NULL,
        affected=NULL, is.deceased=NULL, is.sab.or.top=NULL, is.proband=NULL,
        is.adopted=NULL, are.twins=NULL, are.consanguineous=NULL,
        text.inside.symbol=NULL, text.beside.symbol=NULL,
        text1.below.symbol=NULL, text2.below.symbol=NULL,
        text3.below.symbol=NULL, text4.below.symbol=NULL,
        filename=NULL, main=NULL, device="plot", res=600, ...)

switchPlotfun(method)

}
\arguments{

  \item{family}{
    (Optional) character or numeric vector specifying the family id.
  }

  \item{individual}{
    (Required) character or numeric vector with the ids of the
    individuals.
  }

  \item{father}{
    (Required) character or numeric vector with the id of the father for
    each individual.
  }

  \item{mother}{
    (Required) character or numeric vector with the id of the mother for
    each individual.
  }

  \item{gender}{
    (Required) character, factor or numeric vector specifying the
    gender, with \code{1} or any string starting with \code{"m"} coding
    for male and \code{2} or any string starting with \code{"f"} for
    male; \code{NA} codes for unknown.
  }

  \item{affected}{
    (Optional) numeric or logical vector specifying if the individual is
    affected, \code{0} or \code{FALSE} for not affected, \code{1} or
    \code{TRUE} for affected, \code{NA} for not phenotyped.
  }

  \item{is.deceased}{
    (Optional) numeric of logical vector specifying whether the
    individual is deceased.
  }

  \item{is.sab.or.top}{
    (Optional) numeric or logical vector specifying if the individual is
    the result of a spontaneous abortion or termination of pregnancy.
  }

  \item{is.proband}{
    (Optional) numeric or logical vector specifying whether the
    individual is declared as proband (i.e. the first affected family
    member coming to medical attention).
  }

  \item{is.adopted}{
    (Optional) numeric or logical vector specifying if the individual
    has been adopted.
  }

  \item{are.twins}{
    (Optional) character vector spefifying twins in the
    family. Individuals sharing the same string are recognized as
    twins. The string has to start either with \code{"m_"} or
    \code{"d_"} for monozygotic or dizygotic twins, followed by the
    unique identifier for the twins.
  }

  \item{are.consanguineous}{
    (Optional) character vector specifying cosanguineous couples.
  }

  \item{text.inside.symbol}{
    (Optional) character vector with text to place inside symbols.
  }

  \item{text.beside.symbol}{
    (Optional) character vector with text to place beside symbols.
  }

  \item{text1.below.symbol}{
    (Optional) character vector with text to place below symbols.
  }

  \item{text2.below.symbol}{
    (Optional) character vector with text to place below symbols.
  }

  \item{text3.below.symbol}{
    (Optional) character vector with text to place below symbols.
  }

  \item{text4.below.symbol}{
    (Optional) character vector with text to place below symbols.
  }

  \item{filename}{
    (Optional) file name for the plot. If not specified the result is
    plotted to a temporary file.
  }

  \item{main}{
    (Optional) the title of the plot.
  }

  \item{device}{
    The format of the output file. Can be \code{"ps"}, \code{"pdf"},
    \code{"svg"}, \code{"png"} or \code{"txt"} if \code{HaploPainter} is
    used to create the plot, or \code{"pdf"}, \code{"png"} or
    \code{"plot"} if \code{kinship2} is used for plotting. Note: if
    \code{"plot"} is specified the plot is displayed instead of exported
    to a file.
  }

  \item{res}{
    (Optional) the resolution of the image if saved to a bitmap device.
  }

  \item{method}{
    The method which should be used for plotting, either
    \code{"ks2paint"} (uses \code{kinship2} for plotting) or
    \code{"haplopaint"} (uses \code{HaploPaint}). If not specified, the
    functions switches between the methods.
  }

  \item{...}{
    For \code{plotPed}: additional arguments submitted to the plotting
    function \code{doPlotPed}.
  }


}
\details{
  %% The \code{plotPed} method is simply a wrapper for the \code{doPlotPed}
  %% function, thus all arguments can be used in both. Note that the
  %% \code{plotPed} method for \code{\linkS4class{FAData}} and
  %% \code{\linkS4class{FAResult}} objects are different, i.e. do not
  %% require all the input arguments above, but use the data stored within
  %% the objects to draw the pedigree.

  All arguments passed to the \code{doPlotPed} function have to have the
  same length (with the exception of arguments \code{filename},
  \code{main}, \code{device} and \code{res}) and their order has to
  match the order of the individuals.

  The arguments of the \code{doPlotPed} function represent the input
  parameters supported by \code{HaploPainter}; for more information
  about \code{HaploPainter} refer to http://haplopainter.sourceforge.net/.

  By default, \code{doPlotPed} uses the \code{kinship2} package for
  plotting, but can also be configured by the \code{switchPlotFun} to
  use \code{HaploPainter} instead.
  \code{HaploPainter} is a perl script/tool for pedigree plotting bundled in
  the package that requires however some dependencies that might not be
  present on every system. Thus, the package checks on startup whether
  all requirements for \code{HaploPainter} are available.

  If \code{HaploPainter} is used, the plot can only be exported to a pdf
  or png device, while, if \code{kinship2} is used, the plot can also be
  directly plotted and displayed (if \code{device="plot"} is specified).

  \code{HaploPainter} plotting supports also \code{device = "txt"} in
  which case the pedigree data is exported (in the HaploPainter file
  format) as a tabulator delimited file - no plot is created, the
  name of the file is returned.
  
  Also, the arguments of this function match the arguments for
  \code{HaploPainter} and not all settings can be directly matched to
  settings in \code{kinship2} plotting. The list below lists all
  arguments specific to \code{HaploPainter} and how, if at all, they are
  displayed or used in \code{kinship2} plotting:
  \describe{
    \item{is.sab.or.top}{
      Not supported yet.
    }
    \item{is.proband}{
      The id of individuals which are marked as probands are highlighted
      in red.
    }
    \item{is.adopted}{
      Not supported yet.
    }
    \item{are.twins}{
      Not supported yet.
    }
    \item{are.consanguineous}{
      Not supported yet.
    }
    \item{text2.below.symbol}{
      The text will be plotted on the top left corner of the symbol of
      the respective individual.
    }
    \item{text3.below.symbol}{
      The text will be plotted on the top right corner of the symbol of
      the respective individual.
    }
  }
}
\value{
  \describe{
    \item{switchMethod}{
      A character string representing the plotting function to be used.
    }

    %% \item{plotPed}{
    %%   The name of the file to which the plot was exported.
    %% }

    \item{doPlotPed}{
      The name of the file to which the plot was exported.
    }
  }
}
\author{
  Johannes Rainer.
}
\seealso{
  \code{\link[kinship2]{plot.pedigree}},
  \code{\link{plotPed,FAData-method}},
}
\examples{

## load the Minnesota Breast Cancer record and subset to the
## first families.
data(minnbreast)
family <- minnbreast[minnbreast$famid==4, ]

## draw a pedigree and export it to a pdf file; the file name is
## returned by the function.
doPlotPed(family=family$famid, individual=family$id, father=family$fatherid,
          mother=family$motherid, gender=family$sex, device="pdf")

## switch to the plotting functions of the kinship2 package
switchPlotfun("ks2paint")

## plot the same pedigree, but display it
doPlotPed(family=family$famid, individual=family$id, father=family$fatherid,
          mother=family$motherid, gender=family$sex, device="plot")


}
\keyword{hplot}

