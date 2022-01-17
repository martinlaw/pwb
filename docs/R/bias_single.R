
library(ggplot2)
# Simulate single-stage design:

nsims <- 1e5
N <- 100
theta.vec <- seq(0.1, 0.9, 0.1)
naive <- pwbias <- precision <- numeric(length(theta.vec))
se.list <- vector("list", length(theta.vec))

for (i in 1:length(theta.vec)) {
  responses <- rbinom(nsims, N, theta.vec[i])
  theta.hat <- responses/N
  naive[i] <- mean(theta.hat-theta.vec[i])
  se <- sqrt(theta.hat*(1-theta.hat)/N)
  pwbias[i] <- weighted.mean(x=theta.hat-theta.vec[i], w=1/se^2)
  se.list[[i]] <- se
}

dat <- data.frame(theta=rep(theta.vec, times=2),
                  bias=c(naive, pwbias),
                  type=rep(c("naive", "precision weighted"), each=length(theta.vec)))

ggplot(data=dat, mapping=aes(x=theta, y=bias, col=type))+
  geom_line()+
  labs(title="Bias for single-arm, single-stage trial with N=100")+
    scale_x_continuous(breaks=theta.vec)


# Relative bias:
relative <- data.frame(relative.bias=abs(pwbias)/abs(naive),
                       theta=theta.vec)

ggplot(data=relative, mapping=aes(x=theta, y=relative.bias))+
  geom_line()+
  labs(title="Relative bias (PWB/naive) for single-arm, single-stage trial with N=1000")+
    scale_x_continuous(breaks=theta.vec)+
    geom_hline(aes(yintercept=1), col="grey", linetype="dashed")

# Zoom in:
ggplot(data=relative, mapping=aes(x=theta, y=relative.bias))+
  geom_line()+
  labs(title="Relative bias (PWB/naive) for single-arm, single-stage trial with N=1000")+
  scale_x_continuous(breaks=theta.vec)+
  geom_hline(aes(yintercept=1), col="grey", linetype="dashed")+
  coord_cartesian(ylim=c(-10, 10))

se.df <- data.frame(se=unlist(se.list),
                    theta=rep(as.character(theta.vec), each=nsims))

ggplot(data=se.df, mapping=aes(x=theta, y=se))+
  geom_boxplot()

# Increased variation in SE causes increased variation in weight.
# Increased variation in weight causes weighted mean to be determined by a small
# number of trials with very low SE, which is the trials with theta.hat close to
# 0.5.

hat <- seq(0, 1, 0.01)
seee <- sqrt(hat*(1-hat)/100)
plot(hat, seee)
cbind(hat, seee)
