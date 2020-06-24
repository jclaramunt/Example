#' R package for NETT trees
#'
#' When two treatment alternatives (say A and B) are available  for some problem,
#' one may be interested in qualitative treatment-subgroup interactions. Such
#' interactions imply the existence of subgroups of persons (patients) which are
#' such that in one subgroup Treatment A outperforms Treatment B, whereas the reverse
#' holds in another subgroup. Obviously, this type of interactions is crucial for
#' optimal treatment assignment of future patients. Given baseline characteristics and
#' outcome data from a two-arm Randomized Controlled Trial (RCT), NETT is a tool to identify subgroups that are involved in meaningful
#' qualitative treatment-subgroup interactions. The result of NETT is a tree that
#' partitions the total group of participants (patients) on the basis of their baseline
#' characteristics into three subgroups (i.e., partition classes): Subgroup 1: Those for
#' whom Treatment A is better than Treatment B (P1), Subgroup 2: Those for whom
#' Treatment B is better than Treatment A (P2), and Subgroup 3: Those for whom it does
#' not make any difference (P3).
#'
#' @details \tabular{ll}{
#'   Package: \tab nett\cr
#'   Type: \tab Package\cr
#'   Version: \tab 0.1.0 \cr
#'   Date: \tab 2019-07-05\cr
#'   License: \tab GPL\cr
#' }
#'
#'   The core function of the package is \code{\link{nett}}.
#'
#' @author Maintainer: Elise Dusseldorp <elise.dusseldorp@fsw.leidenuniv.nl>
#' @references Doove, L. L., Dusseldorp, E., Van Deun, K., & Van Mechelen, I. (2015).
#'  A novel method for estimating optimal tree-based treatment regimes in randomized clinical trials.
#'   \emph{Manuscript submitted for publication.}
#' @keywords package
#' @seealso \code{\link{nett}},\code{\link{summary.nett}},\code{\link{nett.control}},
#'   \code{\link{prune.nett}},\code{\link{predict.nett}}
#'
#' @docType package
#' @name nett-package
NULL
