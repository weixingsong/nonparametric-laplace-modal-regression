# From the intermediate results to the entries in Table B.11.

tableB11 = read.csv("InterResult//TableB11//tableB11.csv", skip = 3, header = FALSE)


median=matrix(apply(tableB11,2,median),byrow=T,nrow=12)
iqr=matrix(apply(tableB11,2,IQR),byrow=T,nrow=12)


sink("results/tableB11.txt")

cat("# Simulation Results in Table B.11 \n")
cat("# \n")

cat("# Table B.11 -- median\n")
median
cat("# \n")
cat("# Table B.11 -- interquartile\n")
iqr

sink()