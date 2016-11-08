library(tiger)

symbols=c("BAC","BK","BBT","C","CFG","FITB","HBAN","JPM","KEY","MTB","PNC","STI","USB","WFC","ZION")

data_full=list()

for (i in 1:length(symbols)){
	data=read.csv(paste("C:\\Users\\kushagra\\Documents\\Edu\\Vth Year\\MTP\\snp500data\\snp500data",symbols[i],".csv",sep=""))
	data_ret=diff(data$Close,1);
	data_ret=data_ret-mean(data_ret)
	data_ret=to.uniform(data_ret)
	data_full[[symbols[i]]]=data_ret
	print(mean(data_ret))
}

 write.csv(data_full,"C:\\Users\\kushagra\\Documents\\Edu\\Vth Year\\MTP\\snp500data\\snp500.csv")

library(VineCopula)

d1=data.frame(data_full)
udat=pobs(d1)
pairs(udat)
fit <- RVineStructureSelect(udat)
summary(fit)


