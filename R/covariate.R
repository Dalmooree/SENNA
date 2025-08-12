#' @include create_senna.R
#'
NULL

#' Compute curve parameter
#'
#' Computes curve parameter values (e.g., t) for each spot based on the curve axis (CA) type in the `SENNA` object.
#'
#' @param senna A `SENNA` object
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A `SENNA` object with updated `Spatial` coordinates containing curve parameter
#' @export


GetCurveParam <- function(senna,
                         cores = 1L) {

  if("spline" %in% senna@CurveAxis[["type"]]){
    if("trimmed" %in% senna@CurveAxis[["type"]]){
      if("islet" %in% senna@CurveAxis[["type"]]) {
        pole <- getcvrtisl(senna, cores)
      } else{
        pole <- getcvrtspl(senna, trimmed = TRUE, cores)
      }
    }
    else{
      pole <- getcvrtspl(senna, trimmed = FALSE, cores)
    }

  } else if("straight" %in% senna@CurveAxis[["type"]]) {
    if("trimmed" %in% senna@CurveAxis[["type"]]){
      pole <- getcvrtlin(senna, trimmed = TRUE, cores)
    }
    else{
      pole <- getcvrtlin(senna, trimmed = FALSE, cores)
    }
  }

  senna@Coord[["Spatial"]] <- pole
  return(senna)
}



#' Group curve parameter for downstream analysis
#'
#' Defines the generic function `GroupCurveCovariate()`, used to group spots along the curve parameter for downstream analysis.
#'
#' @importFrom methods setGeneric
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments
#' @return A plot highlighting grouped regions along the curve
#' @export

setGeneric("GroupCurveCovariate", function(senna, ...) {
  standardGeneric("GroupCurveCovariate")
})



#' Grouping curve parameter
#'
#' Groups spatial spots based on their curve parameter (`t`) into discrete bins for downstream analysis.
#'
#' @param senna A `SENNA` object
#' @param bins An integer vector or number of breakpoints for binning
#' @return A `SENNA` object with binned curve parameter added to `Spatial`
#' @export

setMethod("GroupCurveCovariate",
          signature = list(senna = "SENNA"),
          function(senna,
                   bins = NULL){
            if(is.null(bins)) bins <- c(-1e3, senna@CurveAxis[["fun"]][["t"]], 1e3)

            tmp <- dplyr::mutate(
              senna@Coord[["Spatial"]],
              Bin = factor(cut(t, breaks = bins, labels = FALSE, include.lowest = TRUE)))
            tmp <- data.frame(tmp)
            rownames(tmp) <- rownames(senna@Coord[["Spatial"]])

            senna@Coord[["Spatial"]] <- tmp
            return(senna)
          })



#' Get spline curve parameter
#'
#' Computes curve parameter for each spot based on a spline curve axis (CA). Handles both trimmed and untrimmed curves.
#'
#' @param senna A `SENNA` object
#' @param trimmed Logical. Whether the spline is trimmed
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame of updated spatial coordinates with curve parameter values
#' @export

getcvrtspl <- function(senna, trimmed, cores){
  fun <- senna@CurveAxis[["fun"]]
  pole <- senna@Coord[["Spatial"]]

  if(trimmed == TRUE){
    dfun <- trderpolcoef(fun)

    coef1 <- trpol_expan(fun[["x.coef"]])
    coef2 <- trpol_expan(fun[["y.coef"]])
    der1 <- trpol_expan(dfun[["x.deriv"]])
    der2 <- trpol_expan(dfun[["y.deriv"]])
    fun1 <- hfun(coef1, der1, coef2, der2, cores)

    pole <- covt_roots(pole, fun1, coef1, der1, coef2, der2, trimmed = TRUE, cores)
  }

  else{
    dfun <- derpolcoef(fun)

    coef1 <- pol_expan(fun[["x.coef"]])
    coef2 <- pol_expan(fun[["y.coef"]])
    der1 <- pol_expan(dfun[["x.deriv"]])
    der2 <- pol_expan(dfun[["y.deriv"]])
    fun1 <- hfun(coef1, der1, coef2, der2, cores)

    pole <- covt_roots(pole, fun1, coef1, der1, coef2, der2, trimmed = FALSE, cores)
  }

  return(pole)
}



#' Get linear curve parameter
#'
#' Computes curve parameter for each spot based on a linear curve axis (CA). Supports both trimmed and untrimmed modes.
#'
#' @param senna A `SENNA` object
#' @param trimmed Logical. Whether the curve is trimmed
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame of updated spatial coordinates with curve parameter values
#' @export

getcvrtlin <- function(senna, trimmed, cores){
  fun <- senna@CurveAxis[["fun"]]
  pole <- senna@Coord[["Spatial"]]

  if(trimmed == TRUE){
    knot_start <- as.numeric(senna@CurveAxis[["knots"]][1, ])
    knot_end <- as.numeric(senna@CurveAxis[["knots"]][nrow(senna@CurveAxis[["knots"]]), ])

    coef1 <- trpol_expan(fun[["x.coef"]])
    coef2 <- trpol_expan(fun[["y.coef"]])

    pole <- covt_roots_lin(pole, coef1, coef2, trimmed = TRUE,
                           knot_start, knot_end, cores)
  } else{
    coef1 <- pol_expan(fun[["x.coef"]])
    coef2 <- pol_expan(fun[["y.coef"]])

    pole <- covt_roots_lin(pole, coef1, coef2, trimmed = FALSE, cores)
  }
  return(pole)
}



#' Get islet curve parameter
#'
#' Computes curve parameter values for each spot based on a closed islet-type curve axis (CA).
#'
#' @param senna A `SENNA` object
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame of updated spatial coordinates with curve parameter values
#' @export

getcvrtisl <- function(senna, cores){
  fun <- senna@CurveAxis[["fun"]]
  pole <- senna@Coord[["Spatial"]]
  init <- senna@CurveAxis[["knots"]][1,]

  dfun <- trderpolcoef(fun)

  coef1 <- trpol_expan(fun[["x.coef"]])
  coef2 <- trpol_expan(fun[["y.coef"]])
  der1 <- trpol_expan(dfun[["x.deriv"]])
  der2 <- trpol_expan(dfun[["y.deriv"]])
  fun1 <- hfun(coef1, der1, coef2, der2, cores)

  pole <- isl_roots(pole, fun1, coef1, der1, coef2, der2, init, cores)
  return(pole)
}



#' Differentiate curve axis coefficients
#'
#' Computes first derivatives of x and y polynomial coefficients from a spline-based curve axis (CA).
#'
#' @importFrom pracma polyder
#' @param f.list A list containing `t`, `x.coef`, and `y.coef` from a spline curve axis
#' @return A list of differentiated coefficients (`x.deriv`, `y.deriv`) and the corresponding trimmed `t` vector
#' @export

derpolcoef <- function(f.list){
  if(!"t" %in% names(f.list) |
     !"x.coef" %in% names(f.list) |
     !"y.coef" %in% names(f.list)) {
    stop("INVALID SPLINE FUNCTION.")
  }

  t <- f.list$t
  t <- t[2:(length(t)-1)]
  x.c <- f.list$x.coef
  y.c <- f.list$y.coef

  tmp.x <- apply(x.c, 1, pracma::polyder)
  tmp1.x <- lapply(tmp.x, function(l) {
    if(length(l) != 3){
      while(length(l) < 3) {
        l = c(0, l)
      }}
    return(l)
  })
  x.deriv <- t(matrix(unlist(tmp1.x), nrow = 3, byrow = F))

  tmp.y <- apply(y.c, 1, pracma::polyder)
  tmp1.y <- lapply(tmp.y, function(l) {
    if(length(l) != 3){
      while(length(l) < 3) {
        l = c(0, l)
      }}
    return(l)
  })
  y.deriv <- t(matrix(unlist(tmp1.y), nrow = 3, byrow = F))

  return(list(t = t,
              x.deriv = x.deriv,
              y.deriv = y.deriv))
}


#' Differentiate trimmed curve axis coefficients
#'
#' Computes first derivatives of x and y polynomial coefficients from a spline-based trimmed curve axis (CA).
#'
#' @param f.list A list containing `t`, `x.coef`, and `y.coef` from a trimmed spline curve axis
#' @return A list of differentiated coefficients (`x.deriv`, `y.deriv`) and the original `t` vector
#' @export

trderpolcoef <- function(f.list){
  if(!"t" %in% names(f.list) |
     !"x.coef" %in% names(f.list) |
     !"y.coef" %in% names(f.list)) {
    stop("INVALID SPLINE FUNCTION.")
  }

  t <- f.list$t
  x.c <- f.list$x.coef
  y.c <- f.list$y.coef

  tmp.x <- lapply(1:nrow(x.c), function(i){
    der <- pracma::polyder(x.c[i,])
    return(der)
  })
  tmp1.x <- lapply(tmp.x, function(l) {
    if(length(l) != 3){
      while(length(l) < 3) {
        l = c(0, l)
      }}
    return(l)
  })
  x.deriv <- t(matrix(unlist(tmp1.x), nrow = 3, byrow = F))

  tmp.y <- lapply(1:nrow(y.c), function(i){
    der <- pracma::polyder(y.c[i,])
    return(der)
  })
  tmp1.y <- lapply(tmp.y, function(l) {
    if(length(l) != 3){
      while(length(l) < 3) {
        l = c(0, l)
      }}
    return(l)
  })
  y.deriv <- t(matrix(unlist(tmp1.y), nrow = 3, byrow = F))

  return(list(t = t,
              x.deriv = x.deriv,
              y.deriv = y.deriv))
}



#' Expand polynomial coefficients (full curve)
#'
#' Expands segment-wise polynomial coefficients along a full curve axis (CA) for downstream evaluation.
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polypow
#' @importFrom pracma polyadd
#' @param coef_mat A matrix of polynomial coefficients
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A list of expanded polynomial vectors, one per segment
#' @export

pol_expan <- function(coef_mat, cores = 1L) {
  brk <- 1:nrow(coef_mat)
  deg <- 1:ncol(coef_mat)

  l <- parallel::mclapply(brk,
                          function(i) {
                            if(i == 1) {
                              return(coef_mat[1,])
                            }

                            else {
                              vec <- list()
                              for(j in deg){
                                vec[[j]] <- coef_mat[i, j] * pracma::polypow(c(1, -(i - 1)),
                                                                             max(deg)- j)
                              }

                              p_sum <- vec[[1]]

                              for(d in deg[-1]){
                                p_sum <- pracma::polyadd(p_sum, vec[[d]])
                              }

                              return(p_sum)
                            }
                          },
                          mc.cores = cores)

  return(l)
}



#' Expand polynomial coefficients (trimmed curve)
#'
#' Expands segment-wise polynomial coefficients along a trimmed curve axis (CA) for downstream evaluation.
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polypow
#' @importFrom pracma polyadd
#' @param coef_mat A matrix of polynomial coefficients
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A list of expanded polynomial vectors, one per segment
#' @export

trpol_expan <- function(coef_mat, cores = 1L) {
  brk <- 1:nrow(coef_mat)
  deg <- 1:ncol(coef_mat)

  l <- parallel::mclapply(brk,
                          function(i){
                            vec <- list()
                            for(j in deg){
                              vec[[j]] <- coef_mat[i, j] * pracma::polypow(c(1, -i),
                                                                           max(deg)- j)
                            }

                            p_sum <- vec[[1]]

                            for(d in deg[-1]){
                              p_sum <- pracma::polyadd(p_sum, vec[[d]])
                            }

                            return(p_sum)
                          },
                          mc.cores = cores)

  return(l)
}



#' Evaluate h(t) polynomial terms
#'
#' Computes h(t) = x(t)x'(t) + y(t)y'(t) for each segment using polynomial coefficients and their derivatives.
#'
#' @importFrom pracma polyadd
#' @importFrom pracma polymul
#' @importFrom parallel mclapply
#' @param xcoef A list of x polynomial coefficients
#' @param xder A list of x polynomial derivative coefficients
#' @param ycoef A list of y polynomial coefficients
#' @param yder A list of y polynomial derivative coefficients
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A list of h(t) polynomial expressions per segment
#' @export

hfun <- function(xcoef, xder, ycoef, yder, cores = 1L) {
  ind <- 1:max(length(xcoef))
  f <- parallel::mclapply(ind,
                          function(i) {
                            fi <- pracma::polyadd(pracma::polymul(xcoef[[i]], xder[[i]]),
                                                  pracma::polymul(ycoef[[i]], yder[[i]]))
                            return(fi)
                          },
                          mc.cores = cores)

  return(f)
}



#' Find curve parameter for a single point (full spline)
#'
#' Estimates the curve parameter (t) for a spatial point by solving h(t) = ⟨r - γ(t), γ'(t)⟩ = 0 along a full spline curve axis (CA).
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polyadd
#' @importFrom pracma polyval
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics which.min
#' @param point A 2D spatial coordinate (x, y)
#' @param hfun A list of h(t) polynomial terms for each curve segment
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A numeric vector of length 2: estimated curve parameter (t) and squared distance
#' @export

scan_a_cell <- function(point, hfun, xcoef, xder, ycoef, yder, cores = 1L) {
  ind <- 1:max(length(xcoef))

  h <- parallel::mclapply(ind,
                          function(j) {
                            p <- pracma::polyadd(point[[1]] * xder[[j]], point[[2]] * yder[[j]])
                            hi <- pracma::polyadd(hfun[[j]], -p)

                            hrt <- pracma::roots(hi)
                            hrt <- Re(hrt[abs(Im(hrt)) < 1e-5])

                            if(j == 1){
                              hrt <- hrt[hrt < 1 + 1e-5]
                              t <- hrt
                              d <- (point[[1]] - pracma::polyval(xcoef[[j]], hrt)) **2 +
                                (point[[2]] - pracma::polyval(ycoef[[j]], hrt)) **2
                            }
                            else if(j == max(ind)) {
                              hrt <- hrt[hrt >= (j-1 - 1e-5)]
                              t <- hrt
                              d <- (point[[1]] - pracma::polyval(xcoef[[j]], hrt)) **2 +
                                (point[[2]] - pracma::polyval(ycoef[[j]], hrt)) **2
                            }
                            else {
                              hrt <- hrt[hrt >= (j-1 - 1e-5) & hrt < (j + 1e-5)]
                              t <- hrt
                              d <- (point[[1]] - pracma::polyval(xcoef[[j]], hrt)) **2 +
                                (point[[2]] - pracma::polyval(ycoef[[j]], hrt)) **2
                            }

                            all_rts <- BiocGenerics::cbind(t, d)
                            return(all_rts)
                          },
                          mc.cores = cores)

  rts <- do.call(BiocGenerics::rbind, h)
  rts <- rts[BiocGenerics::which.min(rts[ ,2]), ]

  return(rts)
}



#' Find curve parameter for a single point (trimmed spline)
#'
#' Estimates the curve parameter (t) for a spatial point by solving h(t) = ⟨r - γ(t), γ'(t)⟩ = 0 along a trimmed spline curve axis (CA).
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polyadd
#' @importFrom pracma roots
#' @importFrom pracma polyval
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#' @importFrom BiocGenerics which.min
#' @param point A 2D spatial coordinate (x, y)
#' @param hfun A list of h(t) polynomial terms for each curve segment
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A numeric vector of length 2: estimated curve parameter (t) and squared distance. Returns (NaN, NaN) if no solution is found.
#' @export

scan_a_cellt <- function(point, hfun, xcoef, xder, ycoef, yder, cores = 1L) {
  ind <- 1:max(length(xcoef))

  h <- parallel::mclapply(ind,
                          function(j) {
                            p <- pracma::polyadd(point[[1]] * xder[[j]], point[[2]] * yder[[j]])
                            hi <- pracma::polyadd(hfun[[j]], -p)

                            hrt <- pracma::roots(hi)
                            hrt <- Re(hrt[abs(Im(hrt)) < 1e-5])

                            hrt <- hrt[hrt >= (j - 1e-5) & hrt < (j + 1 + 1e-5)]
                            t <- hrt
                            d <- (point[[1]] - pracma::polyval(xcoef[[j]], hrt)) **2 +
                              (point[[2]] - pracma::polyval(ycoef[[j]], hrt)) **2

                            all_rts <- BiocGenerics::cbind(t, d)
                            return(all_rts)
                          },
                          mc.cores = cores)

  rts <- do.call(BiocGenerics::rbind, h)
  if(nrow(rts) == 0){
    rts <- BiocGenerics::rbind(rts, c(NaN, NaN))
  }
  else{
    rts <- rts[BiocGenerics::which.min(rts[ ,2]), ]
  }

  return(rts)
}



#' Find curve parameter for a single point (islet curve)
#'
#' Estimates the curve parameter (t) for a spatial point by solving h(t) = ⟨r - γ(t), γ'(t)⟩ = 0 along an islet-type closed curve axis (CA).  
#' Also considers the distance from the initial knot to enforce closed-curve continuity.
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polyadd
#' @importFrom pracma roots
#' @importFrom pracma polyval
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#' @importFrom BiocGenerics which.min
#' @param point A 2D spatial coordinate (x, y)
#' @param hfun A list of h(t) polynomial terms for each curve segment
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param init The initial knot coordinate used to close
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A numeric vector of length 2: estimated curve parameter (t) and squared distance
#' @export

scan_a_celli <- function(point, hfun, xcoef, xder, ycoef, yder, init, cores = 1L) {
  ind <- 1:max(length(xcoef))

  h <- parallel::mclapply(ind,
                          function(j) {
                            p <- pracma::polyadd(point[[1]] * xder[[j]], point[[2]] * yder[[j]])
                            hi <- pracma::polyadd(hfun[[j]], -p)

                            hrt <- pracma::roots(hi)
                            hrt <- Re(hrt[abs(Im(hrt)) < 1e-5])

                            hrt <- hrt[hrt >= (j - 1e-5) & hrt < (j + 1 + 1e-5)]
                            t <- hrt
                            d <- (point[[1]] - pracma::polyval(xcoef[[j]], hrt)) **2 +
                              (point[[2]] - pracma::polyval(ycoef[[j]], hrt)) **2

                            all_rts <- BiocGenerics::cbind(t, d)
                            return(all_rts)
                          },
                          mc.cores = cores)

  rts <- do.call(BiocGenerics::rbind, h)
  stp <- sum((point[, 1:2] - init[, 1:2])**2)
  rts <- BiocGenerics::rbind(rts, c(1L, stp))
  rts <- rts[BiocGenerics::which.min(rts[ ,2]), ]

  return(rts)
}



#' Find curve parameter for a single point (full linear curve)
#'
#' Estimates the curve parameter (t) for a spatial point by minimizing squared distance to a piecewise linear curve axis (CA), using polynomial root solving.
#'
#' @importFrom parallel mclapply
#' @importFrom pracma polyadd
#' @importFrom pracma polypow
#' @importFrom pracma polyval
#' @importFrom pracma polyder
#' @importFrom pracma roots
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics rbind
#' @importFrom BiocGenerics which.min
#' @param point A 2D spatial coordinate (x, y)
#' @param xcoef Polynomial coefficients for x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A numeric vector of length 2: estimated curve parameter (t) and squared distance
#' @export

scan_a_cell_lin <- function(point, xcoef, ycoef, cores = 1L) {
  ind <- 1:length(xcoef)

  h <- parallel::mclapply(ind,
                          function(j) {
                            if(j == 1){
                              p <- pracma::polyadd(
                                pracma::polypow(pracma::polyadd(xcoef[[j]], c(0, - point[[1]])), 2) ,
                                pracma::polypow(pracma::polyadd(ycoef[[j]], c(0, - point[[2]])), 2)
                              )
                              p <- pracma::polyder(p)
                              rt <- pracma::roots(p)

                              if(rt < j) {
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                              else {
                                rt <- j
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                            }
                            else if(j == max(ind)) {
                              p <- pracma::polyadd(
                                pracma::polypow(pracma::polyadd(xcoef[[j]], c(0, - point[[1]])), 2) ,
                                pracma::polypow(pracma::polyadd(ycoef[[j]], c(0, - point[[2]])), 2)
                              )
                              p <- pracma::polyder(p)
                              rt <- pracma::roots(p)

                              if(rt >= j - 1) {
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                              else {
                                rt <- (j - 1)
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                            }
                            else {
                              p <- pracma::polyadd(
                                pracma::polypow(pracma::polyadd(xcoef[[j]], c(0, - point[[1]])), 2),
                                pracma::polypow(pracma::polyadd(ycoef[[j]], c(0, - point[[2]])), 2)
                              )
                              p <- pracma::polyder(p)
                              rt <- pracma::roots(p)

                              if(rt >= (j - 1 - 1e-5) & rt < (j + 1e-5)) {
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                              else {
                                rt <- ifelse(rt < (j - 1 - 1e-5), j-1, j)
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }
                            }

                            all_rts <- BiocGenerics::cbind(rt, d)
                            return(all_rts)
                          },
                          mc.cores = cores)

  rts <- do.call(BiocGenerics::rbind, h)
  rts <- rts[BiocGenerics::which.min(rts[ ,2]), ]

  return(rts)
}



#' Find curve parameter for a single point (trimmed linear curve)
#'
#' Estimates the curve parameter (t) for a spatial point by minimizing squared distance to a trimmed piecewise linear curve axis (CA), using polynomial root solving.  
#' Returns NaN if the point lies outside the trimmed region defined by `start` and `end`.
#'
#' @importFrom BiocGenerics cbind
#' @importFrom parallel mclapply
#' @importFrom pracma polyadd
#' @importFrom pracma polypow
#' @importFrom pracma polyval
#' @importFrom pracma polyder
#' @importFrom pracma roots
#' @importFrom BiocGenerics rbind
#' @importFrom BiocGenerics which.min
#' @param point A 2D spatial coordinate (x, y)
#' @param xcoef Polynomial coefficients for x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param start Starting coordinate of the trimmed curve
#' @param end Ending coordinate of the trimmed curve
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A numeric vector of length 2: estimated curve parameter (t) and squared distance. Returns (NaN, NaN) if outside trimmed bounds.
#' @export

scan_a_cell_lint <- function(point, xcoef, ycoef, start, end, cores = 1L) {

  if(all((start[1:2] - as.numeric(point[1:2])) %*% c(xcoef[[1]][1], ycoef[[1]][1]) > 0) |
     all((end[1:2] - as.numeric(point[1:2])) %*% c(xcoef[[length(xcoef)]][1], ycoef[[length(ycoef)]][1]) < 0)){

    rt <-
      d <- NaN
    rts <- BiocGenerics::cbind(rt, d)
  }
  else{
    ind <- 1:length(xcoef)
    maxt <- max(ind)

    h <- parallel::mclapply(ind,
                            function(j) {

                              p <- pracma::polyadd(
                                pracma::polypow(pracma::polyadd(xcoef[[j]], c(0, - point[[1]])), 2),
                                pracma::polypow(pracma::polyadd(ycoef[[j]], c(0, - point[[2]])), 2)
                              )
                              p <- pracma::polyder(p)
                              rt <- pracma::roots(p)

                              if(rt >= (j - 1e-5) & rt < (j + 1e-5 + 1)) {
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2
                              }

                              else{
                                rt <- ifelse(rt < (j - 1e-5), j, j + 1)
                                d <- (point[[1]] - pracma::polyval(xcoef[[j]], rt)) **2 +
                                  (point[[2]] - pracma::polyval(ycoef[[j]], rt)) **2

                              }


                              all_rts <- BiocGenerics::cbind(rt, d)
                              return(all_rts)
                            },
                            mc.cores = cores)

    rts <- do.call(BiocGenerics::rbind, h)
    rts <- rts[BiocGenerics::which.min(rts[ ,2]), ]
  }

  return(rts)
}



#' Compute curve parameter and distance for each spot (spline)
#'
#' Computes the curve parameter (t) and Euclidean distance from each spatial spot to the spline-based curve axis (CA).  
#' Supports both full and trimmed spline curves.
#'
#' @importFrom parallel mclapply
#' @importFrom BiocGenerics rbind
#' @importFrom magrittr %>%
#' @param spatial A data frame of spatial coordinates
#' @param hfun A list of h(t) polynomial terms for each curve segment
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param trimmed Logical. Whether the curve is trimmed
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame with additional columns: `t` (curve parameter) and `distance` (Euclidean distance to curve)
#' @export

covt_roots <- function(spatial, hfun, xcoef, xder, ycoef, yder, trimmed, cores = 1L){
  cell_ind <- 1:base::nrow(spatial)

  if(trimmed == TRUE){
    per_cell <- parallel::mclapply(cell_ind, function(i){
      spot <- spatial[i, ]
      percell <- scan_a_cellt(spot, hfun, xcoef, xder, ycoef, yder, cores)
      return(percell)
    },
    mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }
  else if(trimmed == FALSE){
    per_cell <- parallel::mclapply(cell_ind, function(i){
      spot <- spatial[i, ]
      percell <- scan_a_cell(spot, hfun, xcoef, xder, ycoef, yder, cores)
      return(percell)
    },
    mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }


  spatial <- spatial %>%
    mutate(t = res[ ,1], distance = sqrt(res[ ,2]))

  return(spatial)
}



#' Compute curve parameter and distance for each spot (linear curve)
#'
#' Computes the curve parameter (t) and Euclidean distance from each spatial spot to a linear curve axis (CA).  
#' Supports both full and trimmed linear curves.
#'
#' @importFrom parallel mclapply
#' @importFrom BiocGenerics rbind
#' @importFrom magrittr %>%
#' @param spatial A data frame of spatial coordinates
#' @param xcoef Polynomial coefficients for x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param trimmed Logical. Whether the curve is trimmed
#' @param start Starting coordinate of the trimmed curve (required if `trimmed = TRUE`)
#' @param end Ending coordinate of the trimmed curve (required if `trimmed = TRUE`)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame with additional columns: `t` (curve parameter) and `distance` (Euclidean distance to curve)
#' @export

covt_roots_lin <- function(spatial, xcoef, ycoef, trimmed,
                           start = NULL,
                           end = NULL,
                           cores = 1L){
  cell_ind <- 1:nrow(spatial)

  if(trimmed == TRUE){
    per_cell <- parallel::mclapply(cell_ind, function(i){
      spot <- spatial[i, ]
      percell <- scan_a_cell_lint(spot, xcoef, ycoef, start, end, cores)
      return(percell)
    },
    mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }
  else if(trimmed == FALSE){
    per_cell <- parallel::mclapply(cell_ind, function(i){
      spot <- spatial[i, ]
      percell <- scan_a_cell_lin(spot, xcoef, ycoef, cores)
      return(percell)
    },
    mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }


  spatial <- spatial %>%
    mutate(t = res[ ,1], distance = sqrt(res[ ,2]))

  return(spatial)
}



#' Compute curve parameter and distance for each spot (islet curve)
#'
#' Computes the curve parameter (t) and Euclidean distance from each spatial spot to an islet-type closed curve axis (CA), using root solving.
#'
#' @importFrom parallel mclapply
#' @importFrom BiocGenerics rbind
#' @importFrom magrittr %>%
#' @param spatial A data frame of spatial coordinates
#' @param hfun A list of h(t) polynomial terms for each curve segment
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param init The initial knot used to close the curve
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A data frame with additional columns: `t` (curve parameter) and `distance` (Euclidean distance to curve)
#' @export

isl_roots <- function(spatial, hfun, xcoef, xder, ycoef, yder, init, cores = 1L){
  cell_ind <- 1:base::nrow(spatial)

  per_cell <- parallel::mclapply(cell_ind, function(i){
    spot <- spatial[i, ]
    percell <- scan_a_celli(spot, hfun, xcoef, xder, ycoef, yder, init, cores)
    return(percell)
  },
  mc.cores = cores)

  res <- do.call(BiocGenerics::rbind, per_cell)


  spatial <- spatial %>%
    mutate(t = res[ ,1], distance = sqrt(res[ ,2]))

  return(spatial)
}



#' Evaluate spline curve at a given parameter
#'
#' Restores the spatial coordinate by evaluating a cubic spline segment at a given curve parameter `t`, using the provided coefficient matrix.
#'
#' @param t A curve parameter (numeric)
#' @param coef_mat A coefficient matrix for a cubic spline (each row: 4 coefficients for a segment)
#' @return A numeric value representing the spatial coordinate at curve parameter `t`
#' @export

plug_coef <- function(t, coef_mat){
  tmax <- nrow(coef_mat)
  int <- ceiling(t)
  if(int < 1) {
    val <-
      coef_mat[1, 1] * (t - 1)^3 +
      coef_mat[1, 2] * (t - 1)^2 +
      coef_mat[1, 3] * (t - 1) +
      coef_mat[1, 4]
  } else if (int > tmax) {
    val <-
      coef_mat[tmax, 1] * (tmax - 1)^3 +
      coef_mat[tmax, 2] * (tmax - 1)^2 +
      coef_mat[tmax, 3] * (tmax - 1) +
      coef_mat[tmax, 4]
  } else {
    val <-
      coef_mat[int, 1] * (t - (int - 1))^3 +
      coef_mat[int, 2] * (t - (int - 1))^2 +
      coef_mat[int, 3] * (t - (int - 1)) +
      coef_mat[int, 4]
  }

  return(val)
}



#' Evaluate trimmed spline curve at a given parameter
#'
#' Restores the spatial coordinate by evaluating a cubic spline segment at a given curve parameter `t` for a trimmed curve axis.
#'
#' @param t A curve parameter (numeric)
#' @param coef_mat A coefficient matrix for a cubic spline (each row: 4 coefficients for a segment)
#' @return A numeric value representing the spatial coordinate at curve parameter `t`
#' @export

trplug_coef <- function(t, coef_mat){
  int <- floor(t)
  val <-
    coef_mat[int, 1] * (t - int)^3 +
    coef_mat[int, 2] * (t - int)^2 +
    coef_mat[int, 3] * (t - int) +
    coef_mat[int, 4]

  return(val)
}



#' Evaluate linear curve at a given parameter
#'
#' Restores the spatial coordinate by evaluating a piecewise linear curve axis (CA) at a given curve parameter `t`.
#'
#' @param t A curve parameter (numeric)
#' @param coef_mat A coefficient matrix for linear segments (each row: slope and intercept)
#' @return A numeric value representing the spatial coordinate at curve parameter `t`
#' @export

plug_coef_lin <- function(t, coef_mat){
  tmax <- nrow(coef_mat)
  int <- ceiling(t)
  if(int <= 1) {
    val <-
      coef_mat[1, 1] * t +
      coef_mat[1, 2]
  } else if (int >= tmax) {
    val <-
      coef_mat[tmax, 1] * (t - (tmax - 1)) +
      coef_mat[tmax, 2]
  } else {
    val <-
      coef_mat[int, 1] * (t - (int - 1)) +
      coef_mat[int, 2]
  }

  return(val)
}



#' Evaluate trimmed linear curve at a given parameter
#'
#' Restores the spatial coordinate by evaluating a trimmed piecewise linear curve axis (CA) at a given curve parameter `t`.
#'
#' @param t A curve parameter (numeric)
#' @param coef_mat A coefficient matrix for linear segments (each row: slope and intercept)
#' @return A numeric value representing the spatial coordinate at curve parameter `t`
#' @export

trplug_coef_lin <- function(t, coef_mat){
  int <- floor(t)

  val <-
    coef_mat[int, 1] * (t - int) +
    coef_mat[int, 2]

  return(val)
}



