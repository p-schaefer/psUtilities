#' @export

pacleanupFun<-function(data,colname){
  d1<-transpose_fun(data,colname)
  names<-names(which(apply(d1,1,function(x) dim(table(x, useNA = "no"))!=1 && as.numeric(paste0(table(x)[2]))<(sum(as.numeric(paste0(table(x))))-2))))
  d1<-d1[names,]
  return(d1)
}
