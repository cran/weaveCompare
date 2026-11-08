## -----------------------------------------------------------------------------
library(weaveCompare)
a <- outer(0:3, 0:3, function(i,j) as.integer((i-j) %% 4 < 2))
b <- shift_drawdown(a, 2, 1)
c <- a; c[1,1] <- 1-c[1,1]
compare_drawdowns(a,b)[1:7]
duplicate_drawdowns(list(twill=a, shifted=b, edited=c))
d <- drawdown_dist(list(twill=a, shifted=b, edited=c))
as.matrix(d)
plot(hclust(d, method="average"), main="Synthetic structural comparisons")

