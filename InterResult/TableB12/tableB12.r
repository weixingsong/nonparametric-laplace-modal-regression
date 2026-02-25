# From the intermediate results to the entries in Table B.12.


tableB12 = read.csv("InterResult//TableB12//tableB12.csv", skip = 3, header = FALSE)

Bcov=as.matrix(tableB12[1:5])
result=as.matrix(tableB12[6:10])
NBcov=as.matrix(tableB12[11:15])
Nresult=as.matrix(tableB12[16:20])


# MCCL hat_sd and sd

mcclsd=rbind(apply(result,2,sd))
dimnames(mcclsd)=list(c("Stdev"),
                        c("beta0","beta1","beta2","beta3","phi"))

mcclhatsdm=apply(Bcov,2,mean)
mcclhatsds=apply(Bcov,2,sd)


mcclhatsd=rbind(mcclhatsdm,mcclhatsds)
dimnames(mcclhatsd)=list(c("Mean",
                          "Stdev"),
                        c("beta0","beta1","beta2","beta3","phi"))



#Naive Numerical summary of beta and phi

naivesd=rbind(apply(Nresult,2,sd))
dimnames(naivesd)=list(c("Stdev"),
                         c("beta0","beta1","beta2","beta3","phi"))

#Naive Numerical summary of beta and phi  "sd"

naivehatsdm=apply(NBcov,2,mean)
naivehatsds=apply(NBcov,2,sd)

naivehatsd =rbind(naivehatsdm,naivehatsds)
dimnames(naivehatsd)=list(c("Mean",
                           "Stdev"),
                         c("beta0","beta1","beta2","beta3","phi"))

sink("results/tableB12.txt")

cat("# Simulation Results in Table B.12 \n")
cat("# \n")

cat("# Table B.12 -- mcclhatsd\n")
mcclhatsd
cat("# \n")
cat("# Table B.12 -- mcclsd\n")
mcclsd
cat("# \n")
cat("# Table B.12 -- naivehatsd\n")
naivehatsd
cat("# \n")
cat("# Table B.12 -- naivesd\n")
naivesd

sink()

