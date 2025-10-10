filter_matrix_columns <- function(mat, trim_percent = 0.1)
   {
      if (!is.matrix(mat)) mat <- as.matrix(mat)
 
      result <- lapply(1:ncol(mat), function(i) 
         {
            x <- mat[, i]
            q <- quantile(x, probs = c(trim_percent/2, 1 - trim_percent/2), na.rm = TRUE)
            x[x >= q[1] & x <= q[2]]  
          })
  
       max_len <- max(sapply(result, length))
  
      filtered_mat <- sapply(result, function(x)
           {
                length(x) <- max_len
                 x
            })
  
      return(filtered_mat)
   }


     tm0naive <- na.omit(filter_matrix_columns(m0naive, trim_percent = 0.1))
     tm0cor <- na.omit(filter_matrix_columns(m0cor, trim_percent = 0.1))
     tm0huang <- na.omit(filter_matrix_columns(m0huang, trim_percent = 0.1))
     tcorhuang  <- na.omit(filter_matrix_columns(corhuang, trim_percent = 0.1))


    tmod0cormean = apply(tm0cor,2,mean)
    tmod0naivemean = apply(tm0naive,2,mean)
    tmod0huangmean = apply(tm0huang,2,mean)
    tmod0corhuangmean = apply(tcorhuang,2,mean)

    tmisenaivemean = mean(rowMeans((sweep(tm0naive, 2, tmod0naivemean,"-"))^2))
    tmisecormean = mean(rowMeans((sweep(tm0cor, 2, tmod0cormean,"-"))^2))
    tmisehuangmean = mean(rowMeans((sweep(tm0huang, 2, tmod0huangmean))^2))
    tmisecorhuangmean = mean(rowMeans((sweep(tcorhuang, 2, tmod0corhuangmean,"-"))^2))

    tVar <- rbind(tmisenaivemean, tmisecormean,  tmisehuangmean, tmisecorhuangmean)
    rownames(tVar) <- c("Naive", "Corrected", "Huang", "Corhuang")
    tVar


    tmod0cormedian = apply(tm0cor,2,median)
    tmod0naivemedian = apply(tm0naive,2,median)
    tmod0huangmedian = apply(tm0huang,2,median)
    tmod0corhuangmedian = apply(tcorhuang,2,median)

    tmod0cortmean = apply(tm0cor,2,mean,trim=0.05)
    tmod0naivetmean = apply(tm0naive,2,mean,trim=0.05)
    tmod0huangtmean = apply(tm0huang,2,mean,trim=0.05)
    tmod0corhuangtmean = apply(tcorhuang,2,mean,trim=0.05)


     postscript("output.eps", 
           width = 5,   
           height = 5,  
           horizontal = FALSE,  
           paper = "special",   
           onefile = TRUE      
      )

 
   par(mfrow = c(1, 1), mar = c(6, 4, 2, 2) + 0.1, oma = c(1, 0, 0, 0), bty = "l")

   library(latex2exp)    
   x_min <- quantile(xseq,0.1)
   x_max <- quantile(xseq,0.99)
   idx <- xseq >= x_min & xseq <= x_max
   plot(Z,Y,xlab=TeX(r"(log long-term usual intake)"),ylab="log food frequency intake",col = "gray")
   lines(xseq[idx], tmod0naivemedian[idx], type = "l", lwd = 2, col = "green")
   lines(xseq[idx], tmod0huangmedian[idx], lwd = 2, col = "blue")
   lines(xseq[idx], tmod0cormedian[idx], lwd = 2, col = "red")
   lines(xseq[idx], tmod0corhuangmedian[idx], lwd = 2, col = "cyan")

dev.off()

