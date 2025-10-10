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

     corhuang=(1/2)*(m0cor+m0huang)

    tm0naive <- filter_matrix_columns(m0naive, trim_percent = 0.1)
    tm0cor <- filter_matrix_columns(m0cor, trim_percent = 0.1)
    tm0huang <- filter_matrix_columns(m0huang, trim_percent = 0.1)
    tcorhuang  <- filter_matrix_columns(corhuang, trim_percent = 0.1)

    tmisenaivemean = mean(rowMeans((sweep(tm0naive, 2, mod0true,"-"))^2))
    tmisecormean = mean(rowMeans((sweep(tm0cor, 2, mod0true,"-"))^2))
    tmisehuangmean = mean(rowMeans((sweep(tm0huang, 2, mod0true))^2))
    tmisecorhuangmean = mean(rowMeans((sweep(tcorhuang, 2, mod0true,"-"))^2))

    mise <- rbind(tmisenaivemean, tmisecormean,  tmisehuangmean, tmisecorhuangmean)
    rownames(mise) <- c("Naive", "Corrected", "Huang", "Corhuang")
    mise

