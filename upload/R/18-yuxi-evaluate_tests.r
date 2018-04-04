# declare
m0<-0
n_power<-matrix(data=0,nrow=151,ncol=3)
u_power<-matrix(data=0,nrow=151,ncol=3)
c_power<-matrix(data=0,nrow=151,ncol=3)
l_power<-matrix(data=0,nrow=151,ncol=3)

sample<-matrix(NA,nrow=4,ncol=20)
reject<-matrix(data=0,nrow=4,ncol=3)

# cutoff points
alpha<-0.05
t_cut<-qt(alpha,19,lower.tail=FALSE)
b_cut<-qbinom(alpha,20,0.5,lower.tail=FALSE)
w_cut<-qsignrank(alpha,20,lower.tail=FALSE)

# simulate
for (i in 0:150)
{
	m<-3*i/150

	for ( rtrial in 1:1000 )
	{
		# generate random samples
		sample[1,]<-rnorm(20,mean=m,sd=1)
		sample[2,]<-runif(20,min=m-0.5,max=m+0.5)
		sample[3,]<-rcauchy(20,location=m,scale=1)
		for (num in 1:20)
		{
			x<-rexp(1,rate=1)
			y<-rbinom(1,1,0.5)
			if (y)
				x<-x+m
			else
				x<-m-x
			sample[4,num]<-x
		}

		# tests for each distribution
		for (dis in 1:4)
		{
			if ((mean(sample[dis,])-m0)/sqrt(var(sample[dis,])/20)>t_cut)
				reject[dis,1]<-reject[dis,1]+1

			if (sum(sample[dis,]>m0)>b_cut)
				reject[dis,2]<-reject[dis,2]+1

			sample_rank<-rank(abs(sample[dis,]) )
			if ( sum(sample_rank*(sample[dis,]>m0)) > w_cut )
				reject[dis,3]<-reject[dis,3]+1
		}
	}

	# calculate power function
	reject<-reject/1000
	n_power[i+1,]<-reject[1,]
	u_power[i+1,]<-reject[2,]
	c_power[i+1,]<-reject[3,]
	l_power[i+1,]<-reject[4,]

}

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

# Cauchy
plot(xaxis,c_power[,1],type="l",col="red",xlab="theta",ylab="power",ylim=c(0,1))
lines(xaxis,c_power[,2],col="green")
lines(xaxis,c_power[,3],col="blue")
title(main="Cauchy")

