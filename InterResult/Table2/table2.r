# From the intermediate results to the entries in Table 2.


table2 = read.csv("InterResult//Table2//table2.csv", skip = 3, header = FALSE)

Bcov=as.matrix(table2[1:5])
result=as.matrix(table2[6:10])
NBcov=as.matrix(table2[11:15])
Nresult=as.matrix(table2[16:20])


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


sink("results/table2.txt")

cat("# Simulation Results in Table 2 \n")
cat("# \n")

cat("# Table 2 -- mcclhatsd\n")
mcclhatsd
cat("# \n")
cat("# Table 2 -- mcclsd\n")
mcclsd
cat("# \n")
cat("# Table 2 -- naivehatsd\n")
naivehatsd
cat("# \n")
cat("# Table 2 -- naivesd\n")
naivesd

sink()




