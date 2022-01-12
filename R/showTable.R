#' showTable
#'
#' Shows results for single value of theta
#'
#' @param bias.df The data frame containing all results
#' @param theta Response rate
#' @param latex Print Latex code for table (TRUE/FALSE)
#' @param digit Number of rounding digits for bias and SE
#' @return A data frame showing the results for the single theta
#' @import xtable
#' @export
#'
showTable <- function(bias.df, theta, latex=FALSE, digit=4){
  single.theta.df <- bias.df[abs(bias.df$theta-theta)<0.001, ]
  rownames(single.theta.df) <- single.theta.df$type
  single.theta.df$bias <- round(single.theta.df$bias, digit)
  single.theta.df$mean.SE <- round(single.theta.df$mean.SE, digit)
  single.theta.df$emp.SE <- round(single.theta.df$emp.SE, digit)
    if(latex==TRUE){
    print(xtable::xtable(single.theta.df[, 1:4],
                         caption=paste("Bias, mean SE and emp SE when theta=", theta, sep=""),
                         align=rep("r", 5),
                         digits=c(0, 0, digit, digit, digit),
                         label="tab:bias"),
          include.rownames=TRUE,
          booktabs=TRUE,
          sanitize.text.function=function(x){x})
      }
  single.theta.df
}
