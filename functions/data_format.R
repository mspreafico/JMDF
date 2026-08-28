# Function for formatting dataset
formatting.data = function(data){
  data$id <- factor(data$id)
  data$adherent <- factor(data$adherent)
  data$sex <- factor(data$sex)
  data$age <- as.double(data$age)/10
  data$ncom <- as.double(data$ncom)
  return(data)
}


formatting.data.ng = function(data){
  data$id <- as.double(data$id)
  data$adherent <- as.double(data$adherent)
  data$sex <- ifelse(data$sex=='M',1,0)
  setnames(data,'sex','male')
  data$age <- as.double(data$age)/10
  data$ncom <- as.double(data$ncom)
  return(data)
}
