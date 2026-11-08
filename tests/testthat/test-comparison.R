test_that("binary inputs and scalar flags are validated", {
  for (bad in list(NULL, 1, matrix(NA,2), matrix(2,2), matrix(Inf,2), matrix(character(),0,0), matrix('1',1)))
    expect_error(as_drawdown(bad))
  expect_identical(as_drawdown(matrix(TRUE,1)), matrix(1L,1))
  for (v in list(NA, 1, c(TRUE,FALSE))) expect_error(canonical_drawdown(matrix(1,1), reflect=v))
  expect_error(shift_drawdown(matrix(1,1), .5))
  expect_error(shift_drawdown(matrix(1,1), Inf))
  expect_error(drawdown_dist(list()))
  expect_error(drawdown_dist(list(a=matrix(1,1),a=matrix(0,1))))
  expect_error(canonical_drawdown(matrix(1,1), max_work=0))
  expect_error(compare_drawdowns(matrix(0:1,2),matrix(0:1,2),max_work=1))
})

test_that("repeat reduction respects complete periods", {
  a <- matrix(c(1,0,0,1,1,0),2,3)
  tiled <- a[rep(1:2,3),rep(1:3,2)]
  expect_equal(unname(primitive_repeat(tiled)), a)
  expect_equal(dim(primitive_repeat(matrix(1,8,9))), c(1L,1L))
  expect_equal(dim(primitive_repeat(matrix(c(0,1,0),3,1))), c(3L,1L))
  expect_identical(shift_drawdown(a,2,3),as_drawdown(a))
  expect_equal(shift_drawdown(a,-1,-1),shift_drawdown(a,1,2))
})

# Independent reference: enumerate original-domain permutations, without using
# package repeat reduction, shifting, comparison or canonical functions.
oracle <- function(a,b) {
  stopifnot(identical(dim(a),dim(b)))
  nr <- nrow(a); nc <- ncol(a); ans <- 1
  for (r in seq_len(nr)) for (s in seq_len(nc)) {
    rows <- c(seq.int(r,nr),if(r>1) seq_len(r-1) else integer())
    cols <- c(seq.int(s,nc),if(s>1) seq_len(s-1) else integer())
    ans <- min(ans,sum(a != b[rows,cols,drop=FALSE])/length(a))
  }
  ans
}

test_that("all pairs of 2 by 2 matrices match an independent oracle", {
  mats <- lapply(0:15,function(k) matrix(as.integer(intToBits(k))[1:4],2))
  for (a in mats) for (b in mats) {
    expected <- oracle(a,b)
    got <- compare_drawdowns(a,b)
    expect_equal(got$distance,expected)
    expect_identical(got$equivalent,expected==0)
    expect_equal(canonical_drawdown(a)==canonical_drawdown(b),expected==0)
    expect_equal(mean(got$x != got$y),got$distance)
  }
})

test_that("rectangular inputs and optional groups have expected invariants", {
  set.seed(104)
  for (k in seq_len(12)) {
    a <- matrix(sample(0:1,12,TRUE),3,4)
    b <- shift_drawdown(a,1,2)
    expect_equal(compare_drawdowns(a,b)$distance,0)
    expect_equal(canonical_drawdown(a),canonical_drawdown(b))
    expect_equal(compare_drawdowns(a,t(a),transpose=TRUE)$distance,0)
    expect_equal(canonical_drawdown(a,TRUE,TRUE,TRUE),canonical_drawdown(t(1-a),TRUE,TRUE,TRUE))
    expect_equal(compare_drawdowns(a,a[3:1,,drop=FALSE],reflect=TRUE)$distance,0)
    expect_equal(compare_drawdowns(a,1-a,invert=TRUE)$distance,0)
    c <- matrix(sample(0:1,12,TRUE),3,4)
    expect_equal(compare_drawdowns(a,c)$distance,oracle(a,c))
    expect_equal(compare_drawdowns(a,c,TRUE,TRUE,TRUE)$distance,
                 compare_drawdowns(c,a,TRUE,TRUE,TRUE)$distance)
  }
})

test_that("unequal repeats tile onto a least common multiple domain", {
  a <- matrix(c(1,0),2,1); b <- matrix(c(1,0,0),3,1)
  z <- compare_drawdowns(a,b)
  expect_equal(z$domain,c(6,1))
  expect_equal(z$distance,.5)
  expect_equal(compare_drawdowns(matrix(0,1),matrix(1,1))$distance,1)
  expect_false(compare_drawdowns(matrix(0,1),matrix(1,1))$equivalent)
  expect_true(compare_drawdowns(matrix(0,1),matrix(1,1),invert=TRUE)$equivalent)
})

test_that("collection outputs retain labels and exact groups", {
  a <- matrix(c(1,0,0,1),2)
  x <- list(a=a,b=shift_drawdown(a,1),c=matrix(1,1))
  g <- duplicate_drawdowns(x)
  expect_equal(g$group,c(1L,1L,2L))
  expect_equal(g$representative,c('a','a','c'))
  d <- as.matrix(drawdown_dist(x))
  expect_equal(rownames(d),names(x))
  expect_equal(d,t(d)); expect_equal(diag(d),c(a=0,b=0,c=0))
  expect_equal(unname(d[1,]),c(0,0,.5))
  expect_equal(length(drawdown_dist(list(a))),0)
  obj <- structure(list(drawdown=a),class='weave_draft')
  expect_equal(as_drawdown(obj),as_drawdown(a))
})
