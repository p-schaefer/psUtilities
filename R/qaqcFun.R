#' @export

qaqcFun<-function(data,zeros=0.3,maxdiff.zero=0.3,coef=3,type=c("Species","Metrics")) {
  if (!(any(colnames(data)%in%"Site")&any(colnames(data)%in%"Year"))){
    stop("data must contain a column names 'Site' and a column named 'Year'")
  }
  if (length(type)>1){
    stop("Only one type may be selected, either Species or Metrics")
  }

  data$Site<-as.factor(as.character(data$Site))

  if (type=="Species"){
    #zeros is the percent of percent of samples that must be zeros before maxdiff.zero is applied
    #maxdiff.zero is the maximum difference in percent abundance a taxa can change before being flagged as an outlier if there are x zeros
    #coef is how extreme data points are identified. It corresponds to coef*+/-1.58IQR/sqrt(n).
    sites<-unique(data$Site)
    years<-data$Year
    species<-colnames(data)[-c(1:2)]
    flags<-data.frame(matrix(ncol=6))
    colnames(flags)<-c("Class","Site","Year","Species","Value","Median")

    for (i in sites){
      red.data<-data[data$Site==i,]
      rownames(red.data)<-data[data$Site==i,"Year"]

      for (n in species){

        #This block checks for number of zeros and maximum percent abundance differences if greater than the specified number or zeros are present
        if (length(which(red.data[,n]==0))>=zeros*nrow(red.data)){ #check if number of zeros is greater than specified percent
          if ((max(red.data[,n]/rowSums(red.data[,-c(1:2)]))-min(red.data[,n]/rowSums(red.data[,-c(1:2)])))>maxdiff.zero) { #check for any maxdiffs.zero
            yrs<-names(which((red.data[,n]/rowSums(red.data[,-c(1:2)]))>(min(red.data[,n]/rowSums(red.data[,-c(1:2)]))+maxdiff.zero))) #Identify years with maxdiffs.zero
            for (x in yrs){
              entry<-c("maxdiff.zero",paste0(i), x, n, red.data[red.data$Year==x,n], median(red.data[,n]))
              flags<-rbind(flags,entry)
            }
          } else {
            #This is reserved for handling cases where the percent of 0's is greater than the threshold, but the difference in percent abundance is less than the threshold
          }
        }

        #This block looks for outliers in abundance and percent abundance
        if (length(which(red.data[,n]==0))<zeros*nrow(red.data)){
          red.data1<-red.data[,n]
          names(red.data1)<-data[data$Site==i,"Year"]
          pabund<-boxplot.stats(red.data[,n]/rowSums(red.data[,-c(1:2)]),coef=coef)
          raw<-boxplot.stats(red.data1,coef=coef)

          if (!identical(raw$out,numeric(0))){
            yrs<-names(raw$out)
            for (x in yrs){
              entry<-c("abund_outlier",paste0(i), x, n, raw$out[x],raw$stats[3])
              flags<-rbind(flags,entry)
            }
          }

          if (!identical(pabund$out,numeric(0))){
            yrs<-names(pabund$out)
            for (x in yrs){
              entry<-c("percent.abund_outlier",paste0(i), x, n, signif(pabund$out[x],3),signif(pabund$stats[3],3))
              flags<-rbind(flags,entry)
            }
          }
        }
      }


    }

    flags<-flags[-c(1),]
    return(flags)

  }
  if (type == "Metrics") {
    sites<-unique(data$Site)
    years<-data$Year
    metrics<-colnames(data)[-c(1:2)]
    flags<-data.frame(matrix(ncol=6))
    colnames(flags)<-c("Class","Site","Year","Metric","Value","Median")

    for (i in sites){
      red.data<-data[data$Site==i,]
      rownames(red.data)<-data[data$Site==i,"Year"]

      for (n in metrics){
        red.data1<-red.data[,n]
        names(red.data1)<-data[data$Site==i,"Year"]
        raw<-boxplot.stats(red.data1,coef=coef)

        if (!identical(raw$out,numeric(0))){
          yrs<-names(raw$out)
          for (x in yrs){
            entry<-c("metric_outlier",paste0(i), x, n, raw$out[x],raw$stats[3])
            flags<-rbind(flags,entry)
          }
        }

      }
    }
    flags<-flags[-c(1),]
    return(flags)
  }
}
