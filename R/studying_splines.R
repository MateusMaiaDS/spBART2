df <- 5
knots <- NULL
degree = 3
intercept = TRUE
Boundary.knots = range(x)

bs <-
     function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
              Boundary.knots = range(x))
     {
          nx <- names(x)
          x <- as.vector(x)
          nax <- is.na(x)
          if(nas <- any(nax))
               x <- x[!nax]
          if(!missing(Boundary.knots)) {
               Boundary.knots <- sort(Boundary.knots)
               outside <- (ol <- x < Boundary.knots[1]) | (or <- x > Boundary.knots[2])
          } else {
           outside <- FALSE #rep(FALSE, length = length(x))
          }
          ord <- 1 + (degree <- as.integer(degree))
          if(ord <= 1) stop("degree must be integer >= 1")
          if(!missing(df) ) {
               nIknots <- df - ord + (1 - intercept)
               if(nIknots < 0) {
                    nIknots <- 0
                    warning(paste("df was too small; have used ", ord - (1 - intercept)))
                    print("A")
               }
               knots <-
                    if(nIknots > 0) {
                         knots <- seq(from = 0, to = 1, length = nIknots + 2)[-c(1, nIknots + 2)]
                         quantile(x[!outside], knots)
                         print("B")
                    }
          }
          Aknots <- sort(c(rep(Boundary.knots, ord), knots))
          if(any(outside)) {
               warning("Some x values beyond boundary knots may cause ill-conditioned bases")
               derivs <- 0:degree
               scalef <- gamma(1:ord)# factorials
               basis <- array(0, c(length(x), length(Aknots) - degree - 1))
               if(any(ol)) {
                    k.pivot <- Boundary.knots[1]
                    xl <- cbind(1, outer(x[ol] - k.pivot, 1:degree, "^"))
                    tt <- spline.des(Aknots, rep(k.pivot, ord), ord, derivs)$design
                    basis[ol,  ] <- xl %*% (tt/scalef)
               }
               if(any(or)) {
                    k.pivot <- Boundary.knots[2]
                    xr <- cbind(1, outer(x[or] - k.pivot, 1:degree, "^"))
                    tt <- spline.des(Aknots, rep(k.pivot, ord), ord, derivs)$design
                    basis[or,  ] <- xr %*% (tt/scalef)
               }
               if(any(inside <- !outside))
                    basis[inside,  ] <- spline.des(Aknots, x[inside], ord)$design
          }
          else{
               basis <- spline.des(Aknots, x, ord)$design
          }
          if(!intercept){
               basis <- basis[, -1 , drop = FALSE]
          }
          n.col <- ncol(basis)
          if(nas) {
               nmat <- matrix(NA, length(nax), n.col)
               nmat[!nax,  ] <- basis
               basis <- nmat
          }
          dimnames(basis) <- list(nx, 1:n.col)
          a <- list(degree = degree, knots = knots, Boundary.knots =
                         Boundary.knots, intercept = intercept, class = c("bs", "basis"))
          attributes(basis) <- c(attributes(basis), a)
          basis
     }




spline.des <- function(knots, x, ord = 4, derivs = integer(nx))
{
     ## "Design matrix" for a collection of B-splines.  `The' basic function.
     knots <- sort(as.vector(knots))
     x <- as.vector(x)
     nk <- length(knots)
     nx <- length(x)
     ind <- order(x)
     sortx <- x[ind]
     ind <- order(ind)
     if(sortx[1] < knots[ord] || sortx[nx] > knots[nk + 1 - ord]){
          stop(paste("The x data must be in the range", knots[ord], "to",
                     knots[nk + 1 - ord]))
     }
     if(length(derivs)!=nx){
          stop("length of derivs must match length of x")
     }
     ncoef <- nk - ord
     temp <- .C("spline_basis",
                as.double(knots),
                as.integer(ncoef),
                as.integer(ord),
                as.double(sortx),
                as.integer(derivs),
                as.integer(nx),
                design = array(0, c(ord, nx)),
                offsets = integer(nx))
     design <- array(0, c(nx, ncoef))
     d.ind <- array(c(rep(1:nx, rep(ord, nx)),
                      outer(1:ord, temp$offsets, "+")), c(nx * ord, 2))
     design[d.ind] <- temp$design
     list(knots = knots, order = ord, derivs = derivs, design = design[ind,  ])
}

# Putting the arguments
knots <- Aknots
x = x
ord = 4L
derivs = 0L
outer.ok =  FALSE
sparse = FALSE


function (knots, x, ord = 4L, derivs = 0L, outer.ok = FALSE,
          sparse = FALSE)
{
     if ((nk <- length(knots <- as.numeric(knots))) <= 0){
          stop("must have at least 'ord' knots")
     }
     if (is.unsorted(knots)){
          knots <- sort.int(knots)
     }
     x <- as.numeric(x)
     nx <- length(x)
     if (length(derivs) > nx){
          stop("length of 'derivs' is larger than length of 'x'")
     }
     if (length(derivs) < 1L){
          stop("empty 'derivs'")
     }
     ord <- as.integer(ord)
     if (ord > nk || ord < 1){
          stop("'ord' must be positive integer, at most the number of knots")
     }
     if (!outer.ok && nk < 2 * ord - 1){
          stop(gettextf("need at least %s (=%d) knots", "2*ord -1",
                        2 * ord - 1), domain = NA)
     }
     degree <- ord - 1L

     if (need.outer <- any(x < knots[ord] | knots[nk - degree] <
                           x)) {
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
     jj <- c(outer(1L:ord, attr(temp, "Offsets"), `+`))
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

