# code to calculate the bounds
library(lpSolve)
# The inputs of the function are the observable probabilties
# ex. q000 is the probability of being already dead in control group
#     q101 is the probablity of dying given that you have survived up to that point in treatment

calculate_bounds <- function(q000, q100, q111, q001, q101, q110) {
  # Objective function coefficients
  a <- c(0, 0, 0, 0, 0, 1, 1, 0, 0)
  
  # Constraint coefficients
  b <- c(q000 / q100, q001 / q100, q101 / q100, q110 / q100, 1)
  
  # Coefficient matrix
  A <- matrix(c(
    1, 1, 1, -b[1], 0, 0, 0, 0, 0,
    1, 0, 0, 1 - b[2], 1, 0, 0, 0, 0,
    0, 1, 0, -b[3], 0, 1, 1, 0, 0,
    0, 0, 0, -b[4], 1, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 1
  ), nrow = 5, byrow = TRUE)
  
  # Constraint directions
  dir <- rep('==', length(b))
  
  # Lower bound calculation (minimization)
  sol_lb <- lp(direction = 'min', objective.in = a, const.mat = A, const.dir = dir, const.rhs = t(b))
  lb <- sol_lb$objval
  
  # Upper bound calculation (maximization)
  sol_ub <- lp(direction = 'max', objective.in = a, const.mat = A, const.dir = dir, const.rhs = t(b))
  ub <- sol_ub$objval
  
  return(list(lb = lb, ub = ub))
}

# Example usage:

# q000 <- 0.4
# q110 <- 0.2
# q100 <- 1 - q000 - q110
# q001 <- 0.2
# q111 <- 0.7
# q101 <- 1 - q001 - q111
# result <- calculate_bounds(q000, q100, q111, q001, q101, q110)
# print(result$lb)
# print(result$ub)
