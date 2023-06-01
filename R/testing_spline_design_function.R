function (knots, x, ord = 4L, derivs = 0L, outer.ok = FALSE,
          sparse = FALSE)
{
     if ((nk <- length(knots <- as.numeric(knots))) <= 0)
          stop("must have at least 'ord' knots")
     if (is.unsorted(knots))
          knots <- sort.int(knots)
     x <- as.numeric(x)
     nx <- length(x)
     if (length(derivs) > nx) {
          stop("length of 'derivs' is larger than length of 'x'")
          }
     if (length(derivs) < 1L){
          stop("empty 'derivs'")
     }
     ord <- as.integer(ord)
     if (ord > nk || ord < 1) {
          stop("'ord' must be positive integer, at most the number of knots")
     }
     if (!outer.ok && nk < 2 * ord - 1){
          stop(gettextf("need at least %s (=%d) knots", "2*ord -1",
                        2 * ord - 1), domain = NA)
     }
     degree <- ord - 1L
     if (need.outer <- any(x < knots[ord] | knots[nk - degree] <x)) {
          if (outer.ok) {
               in.x <- knots[1L] <= x & x <= knots[nk]
               if ((x.out <- !all(in.x))) {
                    x <- x[in.x]
                    nnx <- length(x)
               }
               dkn <- diff(knots)[(nk - 1L):1]
               knots <- knots[c(rep.int(1L, degree), seq_len(nk),
                                rep.int(nk, max(0L, ord - match(TRUE, dkn > 0))))]
          }
          else stop(gettextf("the 'x' data must be in the range %g to %g unless you set '%s'",
                             knots[ord], knots[nk - degree], "outer.ok = TRUE"),
                    domain = NA)
     }
     temp <- .Call(C_spline_basis, knots, ord, x, derivs)
     ncoef <- nk - ord
     ii <- if (need.outer && x.out) {
          rep.int((1L:nx)[in.x], rep.int(ord, nnx))
     }
     else rep.int(1L:nx, rep.int(ord, nx))
     jj <- c(outer(1L:ord, attr(temp, "Offsets"), "+"))
     if (sparse) {
          if (is.null(tryCatch(loadNamespace("Matrix"), error = function(e) NULL)))
               stop(gettextf("%s needs package 'Matrix' correctly installed",
                             "splineDesign(*, sparse=TRUE)"), domain = NA)
          if (need.outer) {
               jj <- jj - degree - 1L
               ok <- 0 <= jj & jj < ncoef
               methods::as(methods::new("dgTMatrix", i = ii[ok] -
                                             1L, j = jj[ok], x = as.double(temp[ok]), Dim = c(nx,
                                                                                              ncoef)), "CsparseMatrix")
          }
          else methods::as(methods::new("dgTMatrix", i = ii - 1L,
                                        j = jj - 1L, x = as.double(temp), Dim = c(nx, ncoef)),
                           "CsparseMatrix")
     }
     else {
          design <- matrix(double(nx * ncoef), nx, ncoef)
          if (need.outer) {
               jj <- jj - degree
               ok <- 1 <= jj & jj <= ncoef
               design[cbind(ii, jj)[ok, , drop = FALSE]] <- temp[ok]
          }
          else design[cbind(ii, jj)] <- temp
          design
     }
}
