#' Validate a binary weaving drawdown
#'
#' @param x A nonempty logical or numeric binary matrix. A weaveDraft
#'   weave_draft object is also accepted through its drawdown member.
#' @return An integer matrix, with rows representing picks and columns ends.
#' @export
#' @examples
#' as_drawdown(matrix(c(1, 0, 0, 1), 2))
as_drawdown <- function(x) {
  if (inherits(x, "weave_draft")) x <- x$drawdown
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) ||
      !length(x) || anyNA(x) || any(!is.finite(x)) || any(!x %in% c(0, 1)))
    stop("x must be a nonempty, finite binary matrix.", call. = FALSE)
  matrix(as.integer(x), nrow(x), ncol(x), dimnames = dimnames(x))
}

flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
}

#' Reduce a drawdown to its primitive rectangular repeat
#'
#' Only exact periods dividing the observed dimensions are considered.
#' Incomplete repeats are not inferred. Constant matrices reduce to one cell.
#' @param x A binary drawdown; see [as_drawdown()].
#' @return An integer matrix containing the top-left primitive repeat.
#' @export
#' @examples
#' primitive_repeat(matrix(rep(c(1, 0), 8), 4))
primitive_repeat <- function(x) {
  x <- as_drawdown(x)
  period <- function(n, axis) {
    for (k in seq_len(n)) if (n %% k == 0L) {
      ix <- (seq_len(n) - 1L) %% k + 1L
      y <- if (axis == 1L) x[ix, , drop = FALSE] else x[, ix, drop = FALSE]
      if (all(x == y)) return(k)
    }
  }
  x[seq_len(period(nrow(x), 1L)), seq_len(period(ncol(x), 2L)), drop = FALSE]
}

#' Shift a drawdown cyclically
#'
#' Output cell (i,j) takes input cell (i + row_shift,j + col_shift),
#' with indices wrapped around the matrix. Positive shifts move data up/left.
#' @param x A binary drawdown.
#' @param row_shift,col_shift Finite integer offsets; negative offsets allowed.
#' @return An integer matrix of the original dimensions.
#' @export
#' @examples
#' shift_drawdown(matrix(c(1, 0, 0, 0), 2), 1, 0)
shift_drawdown <- function(x, row_shift = 0, col_shift = 0) {
  x <- as_drawdown(x)
  for (v in list(row_shift, col_shift))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != floor(v))
      stop("Offsets must be finite integers.", call. = FALSE)
  x[((seq_len(nrow(x)) - 1 + row_shift %% nrow(x)) %% nrow(x)) + 1,
    ((seq_len(ncol(x)) - 1 + col_shift %% ncol(x)) %% ncol(x)) + 1, drop = FALSE]
}

variants <- function(x, reflect, transpose, invert) {
  flag(reflect, "reflect"); flag(transpose, "transpose"); flag(invert, "invert")
  out <- list(identity = x)
  if (reflect) {
    out$row_reverse <- x[nrow(x):1, , drop = FALSE]
    out$col_reverse <- x[, ncol(x):1, drop = FALSE]
    out$both_reverse <- x[nrow(x):1, ncol(x):1, drop = FALSE]
  }
  if (transpose) out <- c(out, stats::setNames(lapply(out, t), paste0(names(out), "_transpose")))
  if (invert) out <- c(out, stats::setNames(lapply(out, function(y) 1L-y), paste0(names(out), "_invert")))
  out
}

budget <- function(work, max_work) {
  if (!is.numeric(max_work) || length(max_work) != 1L || !is.finite(max_work) || max_work < 1)
    stop("max_work must be a positive finite number.", call. = FALSE)
  if (!is.finite(work) || work > max_work)
    stop("Comparison exceeds max_work; reduce repeat size or raise the limit explicitly.", call. = FALSE)
}

#' Canonical key for a periodic drawdown
#'
#' Keys encode primitive dimensions and row-major bits. The lexicographically
#' smallest representation over allowed transformations and all cyclic shifts
#' is selected. Keys are lossless strings, not cryptographic hashes.
#' @param x A binary drawdown.
#' @param reflect Allow independent reversal of rows and columns.
#' @param transpose Allow axis exchange. This does not invert bits automatically.
#' @param invert Allow exchange of zero and one (face inversion).
#' @param max_work Maximum estimated number of cell comparisons, guarding
#'   exhaustive search. Default 5e7. Memory and runtime still depend on input.
#' @return A single character key. Equivalence is conditional on the options.
#' @export
#' @examples
#' canonical_drawdown(matrix(c(1, 0, 0, 1), 2))
canonical_drawdown <- function(x, reflect = FALSE, transpose = FALSE,
                               invert = FALSE, max_work = 5e7) {
  x <- primitive_repeat(x)
  vv <- variants(x, reflect, transpose, invert)
  budget(length(vv) * length(x)^2, max_work)
  best <- NULL
  for (y in vv) for (i in 0:(nrow(y)-1L)) for (j in 0:(ncol(y)-1L)) {
    z <- shift_drawdown(y, i, j)
    key <- paste0(nrow(z), "x", ncol(z), ":", paste0(as.integer(t(z)), collapse = ""))
    if (is.null(best) || key < best) best <- key
  }
  best
}

#' Align two periodic drawdowns
#'
#' Minimizes the fraction of unequal cells over allowed transformations of y
#' and cyclic shifts. Primitive repeats are tiled to the least common multiple
#' of each axis, so differently sized repeats can be compared without resizing.
#' Ties use the first transformation, then smallest row and column offsets.
#' Both matrices must use the same pick/end and zero/one conventions.
#' @param x,y Binary drawdowns.
#' @param reflect,transpose,invert Transformation flags; see [canonical_drawdown()].
#' @param max_work Cell-comparison budget; see [canonical_drawdown()].
#' @return A list with distance, similarity, equivalent (exact zero distance),
#'   transformation, offsets applied after transformation and repeat reduction,
#'   domain dimensions, and aligned matrices x and y on that domain.
#' @export
#' @examples
#' a <- matrix(c(1, 0, 0, 1), 2)
#' compare_drawdowns(a, shift_drawdown(a, 1))$equivalent
compare_drawdowns <- function(x, y, reflect = FALSE, transpose = FALSE,
                             invert = FALSE, max_work = 5e7) {
  x <- primitive_repeat(x); y <- primitive_repeat(y)
  vv <- variants(y, reflect, transpose, invert)
  gcd <- function(a, b) { while (b != 0) { z <- b; b <- a %% b; a <- z }; a }
  lcm <- function(a, b) a / gcd(a, b) * b
  domains <- lapply(vv, function(z) c(lcm(nrow(x), nrow(z)), lcm(ncol(x), ncol(z))))
  budget(sum(vapply(seq_along(vv), function(k) prod(domains[[k]]) * length(vv[[k]]), numeric(1))), max_work)
  tile <- function(z, d) z[(seq_len(d[1])-1L) %% nrow(z)+1L,
                           (seq_len(d[2])-1L) %% ncol(z)+1L, drop = FALSE]
  best <- Inf; ans <- NULL
  for (k in seq_along(vv)) {
    z <- vv[[k]]; d <- domains[[k]]; a <- tile(x, d)
    for (i in 0:(nrow(z)-1L)) for (j in 0:(ncol(z)-1L)) {
      b <- tile(shift_drawdown(z, i, j), d); score <- mean(a != b)
      if (score < best) {
        best <- score
        ans <- list(distance = score, similarity = 1-score, equivalent = score == 0,
                    transformation = names(vv)[k], row_shift = i, col_shift = j,
                    domain = d, x = a, y = b)
      }
    }
  }
  ans
}

collection <- function(x) {
  if (!is.list(x) || is.data.frame(x) || !length(x))
    stop("drafts must be a nonempty list of drawdowns.", call. = FALSE)
  nm <- names(x)
  if (is.null(nm)) nm <- as.character(seq_along(x))
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("Draft names must be unique and nonempty.", call. = FALSE)
  stats::setNames(lapply(x, as_drawdown), nm)
}

#' Find exact duplicate structures in a collection
#' @param drafts A nonempty list of binary drawdowns, optionally uniquely named.
#' @param reflect,transpose,invert Transformation flags.
#' @param max_work Cell-comparison budget per canonicalization.
#' @return A data frame of draft names, group numbers, first representative and
#'   canonical keys. Groups follow first occurrence, including singleton groups.
#' @export
#' @examples
#' a <- matrix(c(1, 0, 0, 1), 2)
#' duplicate_drawdowns(list(first = a, shifted = shift_drawdown(a, 1)))
duplicate_drawdowns <- function(drafts, reflect = FALSE, transpose = FALSE,
                                invert = FALSE, max_work = 5e7) {
  drafts <- collection(drafts)
  keys <- vapply(drafts, canonical_drawdown, character(1), reflect = reflect,
                 transpose = transpose, invert = invert, max_work = max_work)
  group <- match(keys, unique(keys))
  data.frame(draft = names(drafts), group = unname(group),
             representative = names(drafts)[match(keys, keys)], key = unname(keys),
             row.names = NULL)
}

#' Pairwise periodic Hamming distances
#' @param drafts A nonempty list of binary drawdowns.
#' @param reflect,transpose,invert Transformation flags.
#' @param max_work Cell-comparison budget per pair.
#' @return A stats dist object. Distances ignore phase and optionally other
#'   transformations. Distinct input matrices can have zero distance. No
#'   mechanical, perceptual or cultural meaning is implied. Clustering can use
#'   stats::hclust; Ward clustering is not recommended without Euclidean checks.
#' @export
#' @examples
#' a <- matrix(c(1, 0, 0, 1), 2)
#' drawdown_dist(list(a, matrix(1, 2, 2)))
drawdown_dist <- function(drafts, reflect = FALSE, transpose = FALSE,
                          invert = FALSE, max_work = 5e7) {
  drafts <- collection(drafts)
  flag(reflect, "reflect"); flag(transpose, "transpose"); flag(invert, "invert")
  budget(0, max_work)
  n <- length(drafts); d <- matrix(0, n, n, dimnames = list(names(drafts), names(drafts)))
  if (n > 1L) for (i in seq_len(n-1L)) for (j in (i+1L):n) {
    d[i,j] <- compare_drawdowns(drafts[[i]], drafts[[j]], reflect, transpose, invert, max_work)$distance
    d[j,i] <- d[i,j]
  }
  stats::as.dist(d)
}
