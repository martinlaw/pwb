# showTable: shows results for single value of theta
#
# Arguments:
#
# theta (required)
# latex: show Latex code for table (TRUE/FALSE)
# digit: number of rounding digits for bias and SE
#
showTable <- function(bias.df, theta, latex=FALSE, digit=4){
  single.theta.df <- bias.df[abs(bias.df$theta-theta)<0.001, ]
  rownames(single.theta.df) <- single.theta.df$type
  if(latex==TRUE){
    print(xtable(single.theta.df[, 1:4],
               caption=paste("Bias, mean SE and emp SE when theta=", theta, sep=""),
               align=rep("r", 5),
               digits=c(0, 0, digit, digit, digit),
               label="tab:bias"),
        include.rownames=TRUE,
        booktabs=TRUE,
        sanitize.text.function=function(x){x})
  }
  single.theta.df$bias <- round(single.theta.df$bias, digit)
  single.theta.df$mean.SE <- round(single.theta.df$mean.SE, digit)
  single.theta.df$emp.SE <- round(single.theta.df$emp.SE, digit)
  single.theta.df
}
