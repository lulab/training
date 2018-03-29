# x-axis
xaxis=seq(0,3,by=3/150)

# normal
plot(xaxis,n_power[,1],type="l",col="red",xlab="theta",ylab="power",ylim=c(0,1))
lines(xaxis,n_power[,2],col="green")
lines(xaxis,n_power[,3],col="blue")
title(main="normal")
legend(2,0.3,c("T test", "B test", "W test"),col = c(2,3,4),lty = c(1, 1, 1))

# uniform
plot(xaxis,u_power[,1],type="l",col="red",xlab="theta",ylab="power",ylim=c(0,1))
lines(xaxis,u_power[,2],col="green")
lines(xaxis,u_power[,3],col="blue")
title(main="uniform")
legend(2,0.3,c("T test", "B test", "W test"),col = c(2,3,4),lty = c(1, 1, 1))


# Cauchy
plot(xaxis,c_power[,1],type="l",col="red",xlab="theta",ylab="power",ylim=c(0,1))
lines(xaxis,c_power[,2],col="green")
lines(xaxis,c_power[,3],col="blue")
title(main="Cauchy")
legend(2,0.3,c("T test", "B test", "W test"),col = c(2,3,4),lty = c(1, 1, 1))


# Laplace
plot(xaxis,l_power[,1],type="l",col="red",xlab="theta",ylab="power",ylim=c(0,1))
lines(xaxis,l_power[,2],col="green")
lines(xaxis,l_power[,3],col="blue")
title(main="Laplace")
legend(2,0.3,c("T test", "B test", "W test"),col = c(2,3,4),lty = c(1, 1, 1))
