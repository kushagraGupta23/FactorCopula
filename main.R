model <- mvdc(normalCopula(0.75,dim=4), c("unif","unif","unif","unif"),list(list(min = 0, max =1), list(min = 0, max =1), list(min = 0, max =1), list(min = 0, max =1)))
samples <- rMvdc(1000, mv.NE)

write.csv(t(samples),"C:\\Users\\kushagra\\Documents\\Edu\\Vth Year\\MTP\\kar ke dikhayenge\\datat.csv")

data<-read.csv("C:\\Users\\kushagra\\Documents\\Edu\\Vth Year\\MTP\\kar ke dikhayenge\\datat.csv",row.names=1)

gumbelF<-function(u,v,alpha){
	term1<-(-log(u))^alpha
	term2<-(-log(v))^alpha
	expo<-(term1+term2)^(1/alpha)
	return(exp(-expo))
}

distGumbelF<-function(u,v,alpha){
	eps=.Machine$double.eps^0.50
	umin=0.0
	umax=1.0
	vmin=0.0
	vmax=1.0
	if(u-eps>0){
		umin=u-eps
	}
	if(u+eps<1){
		umax=u+eps
	}
	if(v-eps>0){
		vmin=v-eps
	}
	if(v+eps<1){
		vmax=v+eps
	}
	return((gumbelF(umax,vmax,alpha)-gumbelF(umax,vmin,alpha)-gumbelF(umin,vmax,alpha)+gumbelF(umin,vmin,alpha))/(4*eps*eps))

}

diffGumbelU<-function(u,v,alpha,eps=.Machine$double.eps^0.50){
	return((gumbelF(u+eps,v,alpha)-gumbelF(u-eps,v,alpha))/2*eps)
}

diffGumbelV<-function(u,v,alpha,eps=.Machine$double.eps^0.50){
	return((gumbelF(u,v+eps,alpha)-gumbelF(u,v-eps,alpha))/2*eps)
}

testf<-function(u,v){
	return(gumbelF(u,v,0.71))
}

# Add checks
multply<-function(u,theta,v,FUN){
	mult=1.0
	for( i in 1:length(u)){
		eval=FUN(u[i],v,theta[i])
		if(eval<0){
			print(list(u[i],v,theta[i]))
		}
		mult=mult*(eval);
	}
	return (mult)
}

gaussLegend1D<-function(w,u,v,theta,FUN){
	sum=0.0
	for (i in 1:length(w)){
		sum=sum+w[i]*multply(u,theta,v[i],FUN)
	}
	
	return(sum)
}


# make U list of c()

logLCop<-function(U,theta,FUN){
	l=gaussLegendre(17, 0, 1)
	sum=0.0
	# generate v,w
	v=l$x
	w=l$w
	for(i in 1:length(U)){
		gl=gaussLegend1D(w,t(U[i]),v,theta,FUN)
		#print(gl)
		sum=sum+log(gl);
	}
	return(sum)
}

logLCop(data,c(0.5,0.7,0.9,0.6),distGumbelF)



