     tolrance=0.1
     corhuang=(1/2)*(m0cor+m0huang)




    cvn=(abs(sweep(m0naive, 2, mod0true)) < tolrance*2)
    cvc=(abs(sweep(m0cor, 2, mod0true)) < tolrance*2)
    cvh=(abs(sweep(m0huang, 2, mod0true)) < tolrance*2)
    cvch=(abs(sweep(corhuang, 2, mod0true)) < tolrance*2)


    covercpnaive=mean(apply(cvn,1,all))
    covercpcor=mean(apply(cvc,1,all))
    covercphuang=mean(apply(cvh,1,all))
    covercpch=mean(apply(cvch,1,all))

    coverpnaive=mean(apply(cvn,1,mean))
    coverpcor=mean(apply(cvc,1,mean))
    coverphuang=mean(apply(cvh,1,mean))
    coverpch=mean(apply(cvch,1,mean))

    coverpsdnaive=sd(apply(cvn,1,mean))
    coverpsdcor=sd(apply(cvc,1,mean))
    coverpsdhuang=sd(apply(cvh,1,mean))
    coverpsdch=sd(apply(cvch,1,mean))
   
   
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




    Ccovmod0cormean=mean(abs(mod0cormean-mod0true) < tolrance*2)
    Ccovmod0naivemean=mean(abs(mod0naivemean-mod0true) < tolrance*2)
    Ccovmod0huangmean=mean(abs(mod0huangmean-mod0true) < tolrance*2)
    Ccovmod0corhuangmean=mean(abs(mod0corhuangmean-mod0true) < tolrance*2)

    Ccovmod0cortmean=mean(abs(mod0cortmean-mod0true) < tolrance*2)
    Ccovmod0naivetmean=mean(abs(mod0naivetmean-mod0true) < tolrance*2)
    Ccovmod0huangtmean=mean(abs(mod0huangtmean-mod0true) < tolrance*2)
    Ccovmod0corhuangtmean=mean(abs(mod0corhuangtmean-mod0true) < tolrance*2)

    Ccovmod0cormedian=mean(abs(mod0cormedian-mod0true) < tolrance*2)
    Ccovmod0naivemedian=mean(abs(mod0naivemedian-mod0true) < tolrance*2)
    Ccovmod0huangmedian=mean(abs(mod0huangmedian-mod0true) < tolrance*2)
    Ccovmod0corhuangmedian=mean(abs(mod0corhuangmedian-mod0true) < tolrance*2)


    cover=rbind(c(covercpnaive,covercpcor,covercphuang,covercpch),
                c(coverpnaive,coverpcor,coverphuang,coverpch),
                c(coverpsdnaive,coverpsdcor,coverpsdhuang,coverpsdch),
    c(Ccovmod0naivemean,Ccovmod0cormean,Ccovmod0huangmean,Ccovmod0corhuangmean),
    c(Ccovmod0naivetmean,Ccovmod0cortmean,Ccovmod0huangtmean,Ccovmod0corhuangtmean),
    c(Ccovmod0naivemedian,Ccovmod0cormedian,Ccovmod0huangmedian,Ccovmod0corhuangmedian))
    dimnames(cover)=list(c("covercp","coverp","coversd","Ccovmean","Ccovtmean","Ccovmedian"),
                       c("Naive","Corrected","Huang","Corhuang"))

    
    msecormean = mean((mod0cormean-mod0true)^2)
    msenaivemean = mean((mod0naivemean-mod0true)^2)
    msehuangmean = mean((mod0huangmean-mod0true)^2)
    msecorhuangmean = mean((mod0corhuangmean-mod0true)^2)

    misenaivemean = mean(rowMeans((sweep(m0naive, 2, mod0true,"-"))^2))
    misecormean = mean(rowMeans((sweep(m0cor, 2, mod0true,"-"))^2))
    misehuangmean = mean(rowMeans((sweep(m0huang, 2, mod0true))^2))
    misecorhuangmean = mean(rowMeans((sweep(corhuang, 2, mod0true,"-"))^2))

    mise <- rbind(misenaivemean, misecormean,  misehuangmean, misecorhuangmean)
    rownames(mise) <- c("Naive", "Corrected", "Huang", "Corhuang")

    msecortmean = mean((mod0cortmean-mod0true)^2)
    msenaivetmean = mean((mod0naivetmean-mod0true)^2)
    msehuangtmean = mean((mod0huangtmean-mod0true)^2)
    msecorhuangtmean = mean((mod0corhuangtmean-mod0true)^2)
    
    msecormedian = mean((mod0cormedian-mod0true)^2)
    msenaivemedian = mean((mod0naivemedian-mod0true)^2)
    msehuangmedian = mean((mod0huangmedian-mod0true)^2)
    msecorhuangmedian = mean((mod0corhuangmedian-mod0true)^2)
    
    mse = rbind(c(msenaivemean,msecormean,msehuangmean,msecorhuangmean),
    c(msenaivetmean,msecortmean,msehuangtmean,msecorhuangtmean),
    c(msenaivemedian,msecormedian,msehuangmedian,msecorhuangmedian))
    dimnames(mse) = list(c("Mean","tMean","Median"),
                       c("Naive","Corrected","Huang","Corhuang"))

    mse
    mise
    cover


