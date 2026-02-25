# From the intermediate results to the entries in Table 3.


table3 = read.csv("InterResult//Table3//table3.csv", skip = 2, header = FALSE)

mccl1=as.matrix(table3[1:5])
mccl2=as.matrix(table3[6:10])
naive=as.matrix(table3[11:15])


# MCCL hat_sd and sd

mccl1med=apply(mccl1,2,median)
mccl1iqr=apply(mccl1,2,IQR)

mccl1=rbind(mccl1med,mccl1iqr)

mccl2med=apply(mccl2,2,median)
mccl2iqr=apply(mccl2,2,IQR)
mccl2=rbind(mccl2med,mccl2iqr)

naivemed=apply(naive,2,median)
naiveiqr=apply(naive,2,IQR)
naive=rbind(naivemed,naiveiqr)

dimnames(mccl1)=list(c("Median","IQR"),
                      c("beta0","beta1","beta2","beta3","phi"))
dimnames(mccl2)=list(c("Median","IQR"),
                     c("beta0","beta1","beta2","beta3","phi"))
dimnames(naive)=list(c("Median","IQR"),
                     c("beta0","beta1","beta2","beta3","phi"))



sink("results/table3.txt")

cat("# Simulation Results in Table 3 \n")
cat("# \n")

cat("# Table 3 -- mccl1\n")
mccl1
cat("# \n")
cat("# Table 3 -- mccl2\n")
mccl2
cat("# \n")
cat("# Table 3 -- naive\n")
naive

sink()



