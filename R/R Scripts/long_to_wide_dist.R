long_to_wide_distance <- function (df) {
  
  matNames <- sort(unique(as.character(unlist(df[1:2]))))
  colnames(df)[3] <- "value"
  
  # build matrix of 0s
  myMat <- matrix(0, length(matNames), length(matNames), dimnames = list(matNames, matNames))
  
  # fill in upper triangle
  myMat[as.matrix(df[c(1,2)])] <- df$value
  # fill in the lower triangle
  myMat[as.matrix(df[c(2,1)])] <- df$value

  return(as.dist(myMat))
}