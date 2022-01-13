################################################################################
######                                                                   #######
######             Self-tuning reduced Rank Selection                    #######
######                                                                   #######
################################################################################


#' Self-tuning Rank Selection.
#'
#' \code{STRS} estimates the reduced-rank of the coefficient in multivariate
#' linear regressions. It can also be used to select of the number of factors in
#' factor models.
#'
#' @param data_Y A \eqn{n} by \eqn{m} response data matrix.
#' @param data_X A \eqn{n} by \eqn{p} feature data matrix.
#' @param type String. One of \\{"STRS-DB", "STRS-MC" and "SSTRS"\\}. The
#'   default is "STRS-MC". \itemize{ \item "STRS-DB" and "STRS-MC" are for
#'   general dimensional settings; \item "STRS-DB" uses deterministic
#'   expressions for updating lambda while "STRS-MC" updates lambda by
#'   Monte-Carlo simulations. "STRS-DB" is recommended if "STRS-MC" is
#'   computationally burdensome. \item "SSTRS" is a simpler version of "STRS-DB"
#'   when either \eqn{n >> m} or \eqn{n << m}. }
#' @param rank_X An integer, the specified rank of \code{data_X}. If
#'   unspecified, this is estimated from \code{data_X} as the largest \eqn{k}
#'   such that \deqn{\sigma_k(data_X) \ge rank_tol.}
#' @param self_tune Logical. TRUE if iteratively estimate the rank and
#'   FALSE otherwise. The default is TRUE.
#' @param rep_MC An integer. The number of Monte Carlo simulations. Default is
#'   \eqn{200}.
#' @param rank_tol The tolerence level for determining the rank of \code{data_X}
#'   when \code{rank_X} is NULL. Default is \eqn{1e-4}.
#' @param C A numerical constant for the intial lambda. Default is \eqn{2.01}.
#'
#' @return  The estimated rank.
#'
#' @examples
#'   library(STRS)
#'   est_rank <- STRS(Y, X)
#' @export


STRS <- function(data_Y, data_X, type = "STRS-MC", rank_X = NULL, self_tune = TRUE,
                 rep_MC = 200, rank_tol = 1e-4, C = 2.01) {

  n <- dim(data_Y)[1];  m <- dim(data_Y)[2];  p <- dim(data_X)[2]
  if (is.null(rank_X)) {  # Estimate the rank of X
    q <- max(which(svd(data_X)$d >= rank_tol))
  } else {
    q <- rank_X
  }

  losses <- Compute_total_loss(data_Y, data_X)
  P = data_X %*% MASS::ginv(t(data_X) %*% data_X) %*% t(data_X)

  if (type == "STRS-MC") {
    squarePEs_MC <- MCSquarePEs(q, m, rep_MC)
    resid_MC <- mean(stats::rchisq(rep_MC, (n - q) * m))
    prev_lambda <- C * squarePEs_MC[1]
  } else if (type %in% c("STRS-DB", "SSTRS")) {
    prev_lambda <- C * (sqrt(m) + sqrt(q)) ^ 2
  } else {
    cat("Type must be one of STRS-MC, STRS-DB and SSTRS...")
    stop()
  }

  prev_rank <- Est_rank(losses, prev_lambda, n, q, m)

  if (prev_rank != 0 & self_tune) {
    # iteratively estimate the rank when self_tune is TRUE.
    while(1) {
      curr_lambda <- switch(type,
                            "STRS-DB" = Update_lbd_DB(n, q, m, prev_rank,
                                                      simplified = F),
                            "STRS-MC" = Update_lbd_MC(n, q, m, prev_rank,
                                                      squarePEs_MC, resid_MC),
                            "SSTRS" = Update_lbd_DB(n, q, m, prev_rank,
                                                    simplified = T))
      if (curr_lambda >= prev_lambda) {
        cat("The new lambda is smaller than the old one! Need to modify the
            leading constant...\n")
        stop()
      }
      curr_rank <- Est_rank(losses, curr_lambda, n, q, m, prev_rank)
      if (curr_rank == prev_rank)
        return(curr_rank)
      prev_rank <- curr_rank
      prev_lambda <- curr_lambda
    }
  }
  return(prev_rank)
}


