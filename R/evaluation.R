#' Title
#'
#' @param time 
#' @param progression 
#'
#' @return
#' @export
#'
#' @examples
evaluate.trajectory <- function(time, progression) {
  stime <- sort(time)
  
  # if some timepoints contain multiple cells, randomly order these
  diff <- stime[-1]-stime[-length(stime)]
  min.diff <- min(diff[diff!=0])
  noises <- runif(length(time), 0, .01) * min.diff
  noised.time <- time + noises
  rank <- rank(noised.time)
  
  consistencies <- t(sapply(seq_along(rank), function(i) {
    left <- rank[[i]] > rank
    right <- rank[[i]] < rank
    up <- progression[[i]] > progression
    down <- progression[[i]] < progression
    
    leftup <- sum(left & up)
    leftdown <- sum(left & down)
    rightup <- sum(right & up)
    rightdown <- sum(right & down)
    
    one <- (leftup + rightdown)
    two <- (rightup + leftdown)
    
    consistency.direction1 <- one / (one + two) * 2 - 1
    consistency.direction2 <- two / (one + two) * 2 - 1
    
    c(consistency.direction1, consistency.direction2)
  }))
  
  consistency <- colMeans(consistencies)
  
  max(consistency)
}


#' Title
#'
#' @param space 
#' @param progression 
#' @param k 
#'
#' @return
#' @export
#' 
#' @importFrom class knn
#'
#' @examples
evaluate.space <- function(space, progression, k=5) {
  pred <- sapply(seq_len(nrow(space)), function(i) {
    as.integer(as.character(class::knn(space[-i,,drop=F], space[i,,drop=F], progression[-i], k=k)))
  })
  mean(progression==pred)
}