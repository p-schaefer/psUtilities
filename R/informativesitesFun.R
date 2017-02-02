#' @export

informativesitesFun<-function(data,metric,na.break=0.5,zeros.break=0.6) {
  d1<-transpose_fun(data,metric)
  iqr.range<-apply(d1,1,function(x) abs(IQR(x,na.rm=T)))>0
  na.count<-apply(d1,1,function(x) sum(is.na(x)))<(na.break*ncol(d1))
  zeros.count<-!apply(d1,1,function(x) length(which(x==0))>(zeros.break*ncol(d1)))
  informative.sites<-rownames(d1)[!apply(cbind(iqr.range,na.count,zeros.count
  ),1,function(x) any(x==FALSE))]
  return(informative.sites)
}
