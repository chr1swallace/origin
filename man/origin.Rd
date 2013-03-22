\name{origin}
\alias{origin}
\alias{print.origin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Parent of origin testing
}
\description{
Parent of origin testing using the method proposed by Weinberg, with extension to accomodate more than one affeected offspring in a family.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{origin(genotype, data = sys.frame(sys.parent()), pedigree, id,
              id.father, id.mother, affected, codes = 2, first = FALSE,
              verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
    factor with 3 levels or another object which can be coerced to a
    factor with 3 levels describing genotype.  NB if coerced to a
    factor, it is assumed the levels are ordered so that the
    heterozygote appears in the middle
%%     ~~Describe \code{genotype} here~~
}
  \item{data}{
    data.frame containing genotype (and other variables) if not in
    current environment.
    %%     ~~Describe \code{data} here~~
}
\item{pedigree}{
  name of the variable grouping individuals from the same nuclear family
%%     ~~Describe \code{pedigree} here~~
}
\item{id}{
  name of the variable identifying subjects within pedigrees
%%     ~~Describe \code{id} here~~
}
  \item{id.father}{
    id of subject's father (0 indicates this subject has no father in
    the pedigree)
    %%     ~~Describe \code{id.father} here~~
}
  \item{id.mother}{
    id of subject's mother (0 indicates this subject has no mother in
    the pedigree)
%%     ~~Describe \code{id.mother} here~~
}
\item{affected}{
  whether an individual is affected.  See codes.
%%     ~~Describe \code{affected} here~~
}
\item{codes}{
  value(s) of affected variable which indicates subject is affected.  All
  other values (apart from NA) will be assumed unaffected.
%%     ~~Describe \code{codes} here~~
}
\item{first}{
  if \code{TRUE} only the first affected subject with genotype data will
  be used in each pedigree.
%%     ~~Describe \code{first} here~~
}
\item{verbose}{
  if \code{TRUE} print out (hopefully) useful information as data is
  reformatted for the parent of origin test.
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The return value is a list, with class \code{c("list", "origin")}.
  Its entries are:
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{M3}{Model with imprinting + maternal genotype effects}
 \item{M2}{Model with maternal genotype effects}
 \item{M1}{Model with imprinting effect}
 \item{M0}{Model with no effect (for completeness)}
 %% ...
 A print method is defined, which prints a summary of the models and
 their deviance.  To access the individual models, do:
 \code{models<-origin(snp,data=mydata)
   print(models) ## print comparison
   models$M3 ## print full model
   summary(models$M1) ## print imprinting model with some detail
 }
 }
\references{  
  Weinberg et al. Methods for detection of parent-of-origin effects in
  genetic studies of case-parents triads. \emph{Am J Hum Genet}
  65:229-235. 1999.

  Wallace et al. The imprinted DLK1-MEG3 gene region on chromosome
  14q32.2 alters susceptibility to type 1 diabetes. \emph{Nat
    Genet}. 
  \url{http://dx.doi.org/doi:10.1038/ng.493}
  %% ~put references to the literature/web site here ~
}
\author{
  Chris Wallace
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
