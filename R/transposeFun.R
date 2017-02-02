#' @export

transposeFun<- function(data,metric){
  if (!(any(colnames(data)%in%"Site")&any(colnames(data)%in%"Year"))){
    stop("data must contain a column names 'Site' and a column named 'Year'")
  }

  data$Site<-as.factor(as.character(data$Site))

  sites<-sort(unique(as.character(data$Site)))
  years<-sort(unique(as.numeric(data$Year)))
  output<-data.frame(matrix(nrow=length(sites),ncol=length(years)))
  colnames(output)<-years
  rownames(output)<-sites
  for (i in sites) {
    output[i,which(years%in%data[data$Site%in%i,c(metric,"Year")]$Year)]<-data[data$Site%in%i,metric]
  }

  return(output)
}
