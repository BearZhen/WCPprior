\name{Get_WD}
\alias{Get_WD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Get_WD(flexxible_density, lower_bound, upper_bound, para_name_flexible, para_value_flexible, s, p, args_flexible = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flexxible_density}{
%%     ~~Describe \code{flexxible_density} here~~
}
  \item{lower_bound}{
%%     ~~Describe \code{lower_bound} here~~
}
  \item{upper_bound}{
%%     ~~Describe \code{upper_bound} here~~
}
  \item{para_name_flexible}{
%%     ~~Describe \code{para_name_flexible} here~~
}
  \item{para_value_flexible}{
%%     ~~Describe \code{para_value_flexible} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{args_flexible}{
%%     ~~Describe \code{args_flexible} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (flexxible_density, lower_bound, upper_bound, para_name_flexible, 
    para_value_flexible, s, p, args_flexible = list()) 
{
    string <- deparse(substitute(para_name_flexible))
    args_flexible[[string]] = para_value_flexible
    flexible_density_unipara = function(x) {
        args_flexible[["value"]] = x
        return(do.call(flexible_func, args_flexible))
    }
    integrand = function(x) {
        temp = flexible_density_unipara(x) * sum((x - s) * (x - 
            s))^p
        return(abs(temp))
    }
    integral = integrate(Vectorize(integrand), lower_bound, upper_bound, 
        subdivisions = 2000)$value
    return(integral^(1/p))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
