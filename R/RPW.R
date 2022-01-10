#' randomiseRPW
#'
#' Simulates a randomised play-the-winner trial
#'
#' @param N sample size
#' @param Pa response probabaility for arm A
#' @param Pb response probabaility for arm B
#' @return Vector of results: number of patients and responses in arm A and
#' number of patients and responses in arm B.
#' @export
#'
#'
randomiseRPW = function(N, Pa, Pb){

  P = c(Pa, Pb)

  allocation = response = rep(0,N)

  allocation[1] = rbinom(1,1,0.5)

  response[1] = rbinom(1,1,P[allocation[1]+1])

  urn = ifelse(response[1], allocation[1], 1-allocation[1])

  urn = c(urn, rep(NA,N-1))

  for (i in 2:N){

    ball = sample(urn[1:(i-1)], 1) # chose one ball

    allocation[i] = ball

    response[i] = rbinom(1,1,P[allocation[i]+1])

    urn[i] = ifelse(response[i], ball, 1-ball)

  }

  return(c(sum(allocation==0), sum(response[allocation==0]),
           sum(allocation==1), sum(response[allocation==1])))

}
