amplify <-
function (gamma, lambda) 
{
# This is the amplify function from sensitivitymv version 1.3
    stopifnot(length(gamma) == 1)
    stopifnot(gamma > 1)
    stopifnot(min(lambda) > gamma)
    delta <- (gamma * lambda - 1)/(lambda - gamma)
    names(delta) <- lambda
    delta
}
