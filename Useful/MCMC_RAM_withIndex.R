##==============================================================================
## Version of "MCMC" which passes the index of the current MCMC iteration into
## the log-posterior function provided by the user.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================
## Copyright 2016 Tony Wong, Alexander Bakker
## This file is part of BRICK (Building blocks for Relevant Ice and Climate
## Knowledge). BRICK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## BRICK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

MCMC.withIndex <- function (p, n, init, scale = rep(1, length(init)), adapt = !is.null(acc.rate),
    acc.rate = NULL, gamma = 0.5, list = TRUE, n.start = 0, ...)
{
    if (adapt & !is.numeric(acc.rate))
        stop("Argument \"acc.rate\" is missing!")
    if (gamma < 0.5 | gamma > 1)
        stop("Argument \"gamma\" must be between 0.5 and 1!")
    if (is.numeric(adapt))
        n.adapt <- adapt
    if (adapt == TRUE)
        n.adapt <- Inf
    if (adapt == FALSE)
        n.adapt <- 0
    d <- length(init)
    X <- matrix(NA, ncol = d, nrow = n)
    colnames(X) <- names(init)
    X[1, ] <- init
    p.val <- rep(NA, n)
    p.val[1] <- p(X[1, ], ...)
    if (length(scale) == d) {
        M <- diag(scale)
    }
    else {
        M <- scale
    }
    if (ncol(M) != length(init))
        stop("Length or dimension of 'init' and 'scale' do not match!")
    if (length(init) == 1)
        stop("One-dimensional sampling is not possible!")
    S <- t(chol(M))
    cat("  generate", n, "samples \n")
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    update.step <- max(5, floor(n/100))
    k <- 0
    for (i in 2:n) {
        if (i%%update.step == 0) {
            setTxtProgressBar(pb, i)
        }
        U <- rnorm(d)
        X.prop <- c(X[i - 1, ] + S %*% U)
        names(X.prop) <- names(init)
        p.val.prop <- p(X.prop, index=i, ...)
        alpha <- min(1, exp(p.val.prop - p.val[i - 1]))
        if (!is.finite(alpha))
            alpha <- 0
        if (runif(1) < alpha) {
            X[i, ] <- X.prop
            p.val[i] <- p.val.prop
            k <- k + 1
        }
        else {
            X[i, ] <- X[i - 1, ]
            p.val[i] <- p.val[i - 1]
        }
        ii <- i + n.start
        if (ii < n.adapt) {
            adapt.rate <- min(5, d * ii^(-gamma))
            M <- S %*% (diag(d) + adapt.rate * (alpha - acc.rate) *
                U %*% t(U)/sum(U^2)) %*% t(S)
            eig <- eigen(M, only.values = TRUE)$values
            tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
            if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) >
                tol)) {
                M <- as.matrix(Matrix::nearPD(M)$mat)
            }
            S <- t(chol(M))
        }
    }
    close(pb)
    acceptance.rate <- round(k/(n - 1), 3)
    if (list) {
        return(list(samples = X, log.p = p.val, cov.jump = M,
            n.sample = n, acceptance.rate = acceptance.rate,
            adaption = adapt, sampling.parameters = list(sample.density = p,
                acc.rate = acc.rate, gamma = gamma)))
    }
    else {
        cat("Acceptance rate:", acceptance.rate, "\n")
        return(X)
    }
}
##==============================================================================
## End
##==============================================================================
