#' Compute the loss.
#'
#' Compute the loss \eqn{||Y-(PY)_r||^2} at different rank \eqn{r}.
#'
#' @param Y A \eqn{n} by \eqn{m} matrix.
#' @param X A \eqn{n} by \eqn{p} matrix.
#'
#' @return A list with two elements: \itemize{
#'   \item a vector of all losses,
#'   \item a vector of \eqn{d_j(PY)^2}.
#' }
#' @noRd

Compute_total_loss = function(Y, X) {
  if (is.null(X)) {
    svdPYs = svd(Y, nu=0, nv=0)$d
    cumPYs = cumsum(svdPYs^2)
    null = norm(Y, "F")^2
    allLoss = c(null, null-cumPYs)
  } else {
    PY = X %*% MASS::ginv(crossprod(X)) %*% crossprod(X, Y)
    svdPYs = svd(PY, nu=0, nv=0)$d
    cumPYs = cumsum(svdPYs^2)
    null = norm(PY, "F")^2
    resid = norm(Y - PY, "F")^2
    allLoss = c(null,null-cumPYs) + resid
  }
  return(list(allLoss=allLoss, PYs=svdPYs^2))
}




#' Compute the rank.
#'
#' Compute the rank from the loss of \code{\link{Compute_total_loss}}.
#'
#' @param losses The first element from \code{\link{Compute_total_loss}}.
#' @param lambda The penalty term.
#' @param n The number of rows of \code{data_Y}.
#' @param q The rank of \code{data_X}.
#' @param m The number of columns of \code{data_Y}.
#' @param lowerRank The smallest value of the target rank. The default is \eqn{0}.
#'
#' @return The estimated rank.
#' @noRd


Est_rank <- function(losses, lambda, n, q, m, lowerRank = 0) {
  squarePYs <- losses$PYs
  R <- min(q, m, length(squarePYs))
  allLoss <- losses$allLoss
  counter <- lowerRank
  Rmax <- min(R, trunc((n * m - 1) / lambda))
  while (counter < Rmax) {
    counter <- counter + 1
    tmpLoss <- allLoss[counter] / (n * m - lambda * (counter - 1))
    if (squarePYs[counter] < lambda * tmpLoss)
      return(counter - 1)
  }
  return(Rmax)
}


#' Approximate singualr values.
#'
#' Compute the singular values \eqn{d_j(PE)^2}, for j = 1, ..., \code{q},
#' by Monte Carlo simulations.
#'
#' @param q The rank of \code{data_X}.
#' @param m The number of columns of \code{data_Y}.
#' @param Nsim The number of repetitions in the MC with default equal to \eqn{200}.
#'
#' @return A vector of singular values from 1 to \code{q}.
#' @noRd


MCSquarePEs = function(q, m, Nsim = 200) {
  N <- min(q, m)
  PEs = matrix(0, N, Nsim)
  for (is in 1:Nsim) {
    PEs[ ,is] = svd(matrix(stats::rnorm(q*m), q, m), nu = 0, nv = 0)$d[1:N] ^ 2
  }
  return(apply(PEs, 1, mean))
}


#' Update lambda by MC.
#'
#' Recalculate lambda based on the selected rank. Use Monte Carlo simulations
#' to approximate the value of \eqn{\frac{d_1^2(PE)+...+d_(2r)^2(PE)}{r}} at \code{r}.
#'
#' @param n The number of rows of \code{data_Y}.
#' @param q The rank of \code{data_X}.
#' @param m The number of columns of \code{data_Y}.
#' @param r The selected rank.
#' @param SquarePEs The sequence of singular values of PE computed from \code{\link{MCSquarePEs}}.
#' @param resid The residual value from MC simulations.
#' @param C A numerical constant with default equal to \eqn{1.05}.
#'
#' @return A numerical value of the new lambda.
#' @noRd

Update_lbd_MC = function(n, q, m, r, SquarePEs, resid, C = 1.05) {

  L = min(q,m)
  denom1 = (resid+sum(SquarePEs)-sum(SquarePEs[1:min(2*r,L)]))/SquarePEs[1]+r
  crit1 = C*n*m/denom1
  if (2 * r <= L) {
    if (2 * r <= L - 2) {
      d1 = SquarePEs[2*r+1]
      d2 = SquarePEs[2*r+2]
    } else {
      d1 = d2 = 0
    }
    denom2 = (resid+sum(SquarePEs)-sum(SquarePEs[1:min(2*r,L)]))/(d1+d2)+r
    crit2 = C*n*m / denom2
  } else {
    crit2 = crit1
  }
  return(max(crit1, crit2))
}


#' Update lambda with deterministic bounds.
#'
#' Calculate the new lambda by using deterministic bounds.
#'
#' @param n The number of rows of \code{data_Y}.
#' @param q The rank of \code{data_X}.
#' @param m The number of columns of \code{data_Y}.
#' @param r The selected rank.
#' @param C A numerical constant with default equal to \eqn{0.01}.
#' @param simplified A logical flag indicating whether or not using the simplified
#'   deterministic bounds. The default is FALSE. It is only suitable to set
#'   \code{simplified} to TRUE if n >> m or m << n.
#'
#' @return A numerical value of the new lambda.
#' @noRd

Update_lbd_DB <- function(n, q, m, r, C = 0.01, simplified = F) {

  L <- min(q, m);  H <- max(q, m);  tmp <- (sqrt(q) + sqrt(m))^2

  if (simplified) {
    numer <- n * m
    if (2 * r >= L) {
      denom <- (1 - C) * (n - q) * m / tmp + r
    } else {
      denom <- (1 - C) * (n * m / 2 / H - r) + r
    }

  } else {
    if (2 * r >= L) {
      numer <- n * m * tmp
      denom <- (1 - C) * (n - q) * m + r * tmp
    } else {
      indices1 <- 1 : (2 * r)
      indices2 <- (2 * r + 1) : L
      bound1 = sum((sqrt(H) + sqrt(L - indices1 + 1))^2)
      bound2 = L * H - sum((sqrt(H) - sqrt(indices2 - 1))^2)
      Rt <- min(bound1, bound2)

      Ut <- max(tmp, (sqrt(H) + sqrt(L-2*r))^2 + (sqrt(H) + sqrt(L-2*r-1))^2)
      numer <- n * m * Ut
      denom <- (1 - C) * n * m - Rt + r * Ut
    }
  }
  return(numer / denom)
}

