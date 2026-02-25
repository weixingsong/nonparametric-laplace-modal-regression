# From the intermediate results to the entries in Table 1.

table1 = read.csv("InterResult//Table1//table1.csv", skip = 3, header = FALSE)
median=matrix(apply(table1,2,median),byrow=T,nrow=12)
iqr=matrix(apply(table1,2,IQR),byrow=T,nrow=12)


sink("results/table1.txt")

cat("# Simulation Results in Table 1 \n")
cat("# \n")

cat("# Table 1 -- median\n")
median
cat("# \n")
cat("# Table 1 -- interquartile\n")
iqr

sink()
