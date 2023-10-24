# Function for formatting dataset
formatting.data = function(data){
  data$id <- factor(data$id)
  data$adherent <- factor(data$adherent)
  data$sex <- factor(data$sex)
  data$age <- as.double(data$age)/10
  data$ncom <- as.double(data$ncom)
  return(data)
}