\name{plot.gtobitscad}
\alias{plot.gtobitscad}
\title{Plot coefficients from a gtobitscad object}

\description{
  Plots the coefficient path for different values of \code{lambda} for a fitted \code{gtobitscad} object.
}

\usage{
  \method{plot}{gtobitscad}(x, label = FALSE, ...)
}

\arguments{
  \item{x}{fitted \code{gtobitscad} object.}
  \item{label}{if \code{TRUE}, labels the curves with variable numbers. Default is FALSE.}
  \item{...}{other graphical parameters for \code{plot}}
}

\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{ \code{\link{gtobitscad}} }

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

tscad <- gtobitscad(x = x, y = y)
plot(tscad, label = TRUE)
}
