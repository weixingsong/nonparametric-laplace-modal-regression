    # Yao: ylim=c(0.8,5.0) Huang: ylim=c(-1,6.5)


     corhuang=(1/2)*(m0cor+m0huang)

    mod0cormean = apply(m0cor,2,mean)
    mod0naivemean = apply(m0naive,2,mean)
    mod0huangmean = apply(m0huang,2,mean)
    mod0corhuangmean = apply(corhuang,2,mean)

    mod0cortmean = apply(m0cor,2,mean,trim=0.05)
    mod0naivetmean = apply(m0naive,2,mean,trim=0.05)
    mod0huangtmean = apply(m0huang,2,mean,trim=0.05)
    mod0corhuangtmean = apply(corhuang,2,mean,trim=0.05)

    mod0cormedian = apply(m0cor,2,median)
    mod0naivemedian = apply(m0naive,2,median)
    mod0huangmedian = apply(m0huang,2,median)
    mod0corhuangmedian = apply(corhuang,2,median)



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
  
      filtered_mat <- sapply(result, function(x) {
                                                                            length(x) <- max_len
                                                                              x
                                                                         })
  
       return(filtered_mat)
    }

    tm0naive <- filter_matrix_columns(m0naive, trim_percent = 0.1)
    tm0cor <- filter_matrix_columns(m0cor, trim_percent = 0.1)
    tm0huang <- filter_matrix_columns(m0huang, trim_percent = 0.1)
    tcorhuang  <- filter_matrix_columns(corhuang, trim_percent = 0.1)


    tmod0cormean = apply(tm0cor,2,mean)
    tmod0naivemean = apply(tm0naive,2,mean)
    tmod0huangmean = apply(tm0huang,2,mean)
    tmod0corhuangmean = apply(tcorhuang,2,mean)

    postscript("output.eps", 
           width = 8,   
           height = 8, 
           horizontal = FALSE,  
           paper = "special",   
           onefile = TRUE      
      )
par(mfrow = c(2, 2), mar = c(6, 4, 2, 2) + 0.1, oma = c(1, 0, 0, 0), bty = "l")

    plot(xseq,mod0cormean,type="l",lwd=2,col="red",ylim=c(-0.9,5.0),xlab="x",ylab="mean")
    for(jj in 1:nrow(tm0cor))
     {
      lines(xseq,tm0cor[jj,],col="gray")
     }
    lines(xseq, mod0true,  lwd=2, col="black")
    lines(xseq, mod0naivemean,lwd=2,col="green")
    lines(xseq, mod0cormean, lwd=2, col="red")
    lines(xseq, mod0huangmean,lwd=2,col="blue")
    lines(xseq, mod0corhuangmean,lwd=2,col="cyan")

 legend(x = 0.4,  
        y = 2.5,         
       legend = c("True", "LocPoly", "LocLinear", "LocPL", "Naive"),
       col = c("black", "red", "blue", "cyan", "green"),
       lwd = c(2, 2, 2, 2, 2),
       lty = c(1, 1, 1, 1, 1),
       horiz = FALSE,
       xpd = NA,
       cex = 1.0,  # 缩小字体
       y.intersp = 1.0,   # 缩小文字与线条的垂直间距（默认1.0）
       text.width = NULL,
       bty = "n"
       )

    plot(xseq,mod0cormean,type="l",lwd=2,col="red",ylim=c(-0.9,5.0),xlab="x",ylab="mean")
    for(jj in 1:nrow(tm0cor))
    {
      lines(xseq,tm0huang[jj,],col="gray")
    }
    lines(xseq, mod0true,  lwd=2, col="black")
    lines(xseq, mod0naivemean,lwd=2,col="green")
    lines(xseq, mod0cormean, lwd=2, col="red")
    lines(xseq,mod0huangmean,lwd=2,col="blue")
    lines(xseq,mod0corhuangmean,lwd=2,col="cyan")
 


        
    plot(xseq,mod0cortmean,type="l",lwd=2,col="red",ylim=c(0.8,4.6),xlab="x",ylab="trimmed mean")
    lines(xseq, mod0true,  lwd=2, col="black")
    lines(xseq, mod0naivetmean,lwd=2,col="green")
    lines(xseq, mod0cortmean, lwd=2, col="red")
    lines(xseq,mod0huangtmean,lwd=2,col="blue")
    lines(xseq,mod0corhuangtmean,lwd=2,col="cyan")
       
    
    plot(xseq,mod0cormedian,type="l",lwd=2,col="red",ylim=c(0.8,4.6),xlab="x",ylab="median")
    lines(xseq, mod0true,  lwd=2, col="black")
    lines(xseq, mod0naivemedian,lwd=2,col="green")
    lines(xseq, mod0cormedian, lwd=2, col="red")
    lines(xseq, mod0huangmedian,lwd=2,col="blue")
    lines(xseq, mod0corhuangmedian,lwd=2,col="cyan")
    
dev.off()




