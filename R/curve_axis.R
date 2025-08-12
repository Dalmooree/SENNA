#' @include create_senna.R
#'
NULL


#' Create a curve axis (CA) by connecting the given knots with extrapolation
#'
#' Defines the generic function `FullCurve()`. It is suitable for both 'regionation' scenarios and 'progression' scenarios where the entire FOV is used.
#'
#' @importFrom methods setGeneric
#' @param senna A `SENNA` object
#' @param ... Additional arguments
#' @export

setGeneric("FullCurve", function(senna, ...) {
  standardGeneric("FullCurve")
})



#' Create a Curve Axis (CA) by Connecting the Given Knots without Extrapolation
#'
#' Defines the generic of `TrimmedCurve()`. It is suitable for both 'islet' scenarios and 'progression' scenarios where only part of the FOV is used.
#'
#' @importFrom methods setGeneric
#' @param senna A SENNA object
#' @param ... Additional arguments
#' @export
#'
setGeneric("TrimmedCurve", function(senna, ...) {
  standardGeneric("TrimmedCurve")
})



#' Create a curve axis (CA) by connecting the given knots with extrapolation
#'
#' Defines the method for `FullCurve()` when applied to a `SENNA` object. This function constructs a full curve axis by connecting the input knots and optionally extends the curve beyond the first and last knots. This is suitable for full-field 'progression' or 'regionation' scenarios.
#'
#' @importFrom methods setMethod
#' @param senna A `SENNA` object
#' @param knot_df A data frame of 2D coordinates representing the curve knots
#' @param type A string indicating the curve type; one of `"spline"` or `"straight"`
#' @return A modified `SENNA` object with a fitted full-length curve
#' @export

setMethod("FullCurve",
          signature = list(senna = "SENNA"),
          function(senna,
                   knot_df,
                   type = "spline"){

            t <- seq_along(knot_df$X1)

            if(type == "spline") {
              res <- ctype_spline(senna = senna,
                                  knot_df = knot_df,
                                  t = t)
            } else if(type == "straight") {
              res <- ctype_linear(senna = senna,
                                  knot_df = knot_df,
                                  t = t)
            } else{
              stop("The type must be either 'spline' or 'straight'.")
            }

            return(res)
          })



#' Generate trimmed curve axis
#'
#' Method for the generic `TrimmedCurve()` applied to a `SENNA` object.
#'
#' @importFrom methods setMethod
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @param type Type of curve. Must be one of `'spline'`, `'straight'`, or `'islet'`
#' @return A `SENNA` object with curve axis added
#' @exportMethod TrimmedCurve

setMethod("TrimmedCurve",
          signature = list(senna = "SENNA"),
          function(senna,
                   knot_df,
                   type = "spline"){
            if(type == "spline") {
              t <- seq_along(knot_df$X1)
              res <- trtype_spline(senna = senna,
                                   knot_df = knot_df,
                                   t = t)
            } else if(type == "straight") {
              t <- seq_along(knot_df$X1)
              res <- trtype_linear(senna = senna,
                                   knot_df = knot_df,
                                   t = t)
            } else if(type == "islet"){
              res <- trtype_islet(senna = senna,
                                  knot_df = knot_df)
            } else{
              stop("In `TrimmedCurve()`, the type must be one of 'spline', 'straight', or 'islet'.")
            }

            return(res)
          })



#' Activate specific clusters for analysis
#'
#' Reflects cluster label information in the analysis based on the saved annotation.
#'
#' @param senna A `SENNA` object
#' @param id A vector of cluster identities to activate
#' @return A modified `SENNA` object with activated cluster information
#' @export

ActiveIdent <- function(senna,
                        id) {
  activ <- senna@Gene[["Reference"]][["Annotation"]] %in% id
  sname <- names(senna@Gene[["Reference"]][["Annotation"]][activ])

  senna@Coord[["Activated"]] <- list(ind = TRUE,
                                     ID = id,
                                     cells = sname)

  return(senna)
}



#' Full curve (spline type)
#'
#' Computes spline-based curve axis (CA) coefficients from given knots.
#'
#' @importFrom splines interpSpline
#' @importFrom BiocGenerics rbind
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @param t A numeric vector specifying knot positions
#' @return A `SENNA` object with spline curve axis stored
#' @export

ctype_spline <- function(senna,
                         knot_df,
                         t){
  if(length(t) <= 1){
    stop("To generate curve, at least two knots is needed.")
  }

  else if(length(t) == 2){
    stop("Use straight type for 2 knots.")
  }

  else if(length(t) == 3){
    fx <- spl3k(knot_df$X1)
    fy <- spl3k(knot_df$X2)
    ord <- 2

    t_start <- c(-fx[1,4] / fx[1,3],
                 (1-fx[1,4]) / fx[1,3],
                 -fy[1,4] / fy[1,3],
                 (1-fy[1,4]) / fy[1,3])
    t_start <- max(t_start[t_start <= 1])

    t_end <- c(-fx[nrow(fx),4] / fx[nrow(fx),3] + max(t),
               (1-fx[nrow(fx),4]) / fx[nrow(fx),3] + max(t),
               -fy[nrow(fy),4] / fy[nrow(fy),3] + max(t),
               (1-fy[nrow(fy),4]) / fy[nrow(fy),3] + max(t))
    t_end <- min(t_end[t_end >= max(t)])

    t <- c(t_start, t, t_end)
  }

  else{
    ord <- 3

    fx <- splines::interpSpline(knot_df$X1 ~ t)$coefficients
    Cx <- c(fx[1,1] - fx[1,2], fx[1,2], 0, 0)
    fx <- BiocGenerics::rbind(Cx, fx)
    rownames(fx) <- NULL

    fx <- t(apply(fx, 1, rev))

    fy <- splines::interpSpline(knot_df$X2 ~ t)$coefficients
    Cy <- c(fy[1,1] - fy[1,2], fy[1,2], 0, 0)
    fy <- BiocGenerics::rbind(Cy, fy)
    rownames(fy) <- NULL

    fy <- t(apply(fy, 1, rev))

    t_start <- c(-fx[1,4] / fx[1,3],
                 (1-fx[1,4]) / fx[1,3],
                 -fy[1,4] / fy[1,3],
                 (1-fy[1,4]) / fy[1,3])
    t_start <- ifelse(Inf %in% t_start | -Inf %in% t_start,
                      1,
                      max(t_start[t_start <= 1]))

    t_end <- c(-fx[nrow(fx),4] / fx[nrow(fx),3] + max(t),
               (1-fx[nrow(fx),4]) / fx[nrow(fx),3] + max(t),
               -fy[nrow(fy),4] / fy[nrow(fy),3] + max(t),
               (1-fy[nrow(fy),4]) / fy[nrow(fy),3] + max(t))
    t_end <- ifelse(Inf %in% t_end | -Inf %in% t_end,
                    max(t),
                    min(t_end[t_end >= max(t)]))

    t <- c(t_start, t, t_end)

  }

  fun <- list(t = t, x.coef = fx, y.coef = fy)

  t.list <- list(
    type = "spline",
    knots = knot_df,
    fun = fun,
    ord = ord)

  senna@CurveAxis = t.list

  return(senna)
}



#' Full curve (spline, three-point case)
#'
#' Computes spline coefficients when exactly three knots are provided.
#'
#' @param x A numeric vector of length 3 (knot positions)
#' @return A 4×4 coefficient matrix
#' @export

spl3k <- function(x){
  b1 <- -x[3] / 2 + 2 * x[2] - 3 * x[1] / 2
  coef <- matrix(c(0, 0, b1, x[1] - b1,
                   0, x[2] - x[1] - b1, b1, x[1],
                   0, x[3] - 3*x[2] + 2*x[1] + b1, 2*(x[2]-x[1]) - b1, x[2],
                   0, 0, 2*(x[3] - 2*x[2] + x[1]) + b1, x[3]),
                 nrow = 4,
                 ncol = 4,
                 byrow = TRUE)
  return(coef)
}



#' Full curve (piecewise linear)
#'
#' Computes piecewise linear curve axis (CA) coefficients from given knots.
#'
#' @importFrom BiocGenerics rbind
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @param t A numeric vector specifying knot positions
#' @return A `SENNA` object with linear curve axis stored
#' @export

ctype_linear <- function(senna,
                         knot_df,
                         t){
  ord <- 1
  
  x1 <- knot_df$X1
  fx <- lapply(c(0, t), function(i){
    if(i == 0) {
      vec <- c(x1[i + 2] - x1[i + 1], 2 * x1[i + 1] - x1[i + 2])
    } else if (i == max(t)) {
      vec <- c(x1[i] - x1[i - 1], x1[i])
    } else{
      vec <- c(x1[i + 1] - x1[i], x1[i])
    }
    return(vec)
  })
  fx <- do.call(BiocGenerics::rbind, fx)
  
  x2 <- knot_df$X2
  fy <- lapply(c(0, t), function(i){
    if(i == 0) {
      vec <- c(x2[i + 2] - x2[i + 1], 2 * x2[i + 1] - x2[i + 2])
    } else if (i == max(t)) {
      vec <- c(x2[i] - x2[i - 1], x2[i])
    } else{
      vec <- c(x2[i + 1] - x2[i], x2[i])
    }
    return(vec)
  })
  fy <- do.call(BiocGenerics::rbind, fy)
  
  t_start <- c(-fx[1,2] / fx[1,1],
               (1 - fx[1,2]) / fx[1,1],
               -fy[1,2] / fy[1,1],
               (1 - fy[1,2]) / fy[1,1])
  t_start <- ifelse(Inf %in% t_start | -Inf %in% t_start,
                    1,
                    max(t_start[t_start <= 1]))
  
  t_end <- c(-fx[nrow(fx),2] / fx[nrow(fx),1] + max(t),
             (1-fx[nrow(fx),2]) / fx[nrow(fx),1] + max(t),
             -fy[nrow(fy),2] / fy[nrow(fy),1] + max(t),
             (1-fy[nrow(fy),2]) / fy[nrow(fy),1] + max(t))
  t_end <- ifelse(Inf %in% t_end | -Inf %in% t_end,
                  max(t),
                  min(t_end[t_end >= max(t)]))
  
  t <- c(t_start, t, t_end)
  
  fun <- list(t = t, x.coef = fx, y.coef = fy)
  
  t.list <- list(
    type = "straight",
    knots = knot_df,
    fun = fun,
    ord = ord)
  
  senna@CurveAxis = t.list
  
  return(senna)
}



#' Trimmed curve (spline)
#'
#' Computes spline-based trimmed curve axis (CA) coefficients from given knots.
#'
#' @importFrom splines interpSpline
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @param t A numeric vector specifying knot positions
#' @return A `SENNA` object with trimmed spline curve axis stored
#' @export

trtype_spline <- function(senna,
                          knot_df,
                          t){
  if(length(t) <= 1){
    stop("To generate curve, at least two knots is needed.")
  }

  else if(length(t) == 2){
    stop("Use 'straight' type for 2 knots.")
  }

  else if(length(t) == 3){
    fx <- spl3k(knot_df$X1)
    fx <- fx[2:3,]
    fy <- spl3k(knot_df$X2)
    fy <- fy[2:3,]
    ord <- 2
  }

  else{
    ord <- 3

    fx <- splines::interpSpline(knot_df$X1 ~ t)$coefficients
    fx <- fx[-nrow(fx),]

    fx <- t(apply(fx, 1, rev))

    fy <- splines::interpSpline(knot_df$X2 ~ t)$coefficients
    fy <- fy[-nrow(fy),]

    fy <- t(apply(fy, 1, rev))
  }

  fun <- list(t = t, x.coef = fx, y.coef = fy)

  t.list <- list(
    type = c("spline", "trimmed"),
    knots = knot_df,
    fun = fun,
    ord = ord)

  senna@CurveAxis = t.list

  return(senna)
}



#' Trimmed curve (piecewise linear)
#'
#' Computes piecewise linear trimmed curve axis (CA) coefficients from given knots.
#'
#' @importFrom BiocGenerics rbind
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @param t A numeric vector specifying knot positions
#' @return A `SENNA` object with trimmed linear curve axis stored
#' @export

trtype_linear<- function(senna,
                         knot_df,
                         t){
  ord <- 1

  x1 <- knot_df$X1
  fx <- lapply(t[-length(t)], function(i){
    vec <- c(x1[i + 1] - x1[i], x1[i])
    return(vec)
  })
  fx <- do.call(BiocGenerics::rbind, fx)

  x2 <- knot_df$X2
  fy <- lapply(t[-length(t)], function(i){
    vec <- c(x2[i + 1] - x2[i], x2[i])
    return(vec)
  })
  fy <- do.call(BiocGenerics::rbind, fy)

  fun <- list(t = t, x.coef = fx, y.coef = fy)

  t.list <- list(
    type = c("straight", "trimmed"),
    knots = knot_df,
    fun = fun,
    ord = 1)

  senna@CurveAxis = t.list

  return(senna)
}



#' Trimmed curve (islet)
#'
#' Computes coefficients for a closed curve axis (CA) based on islet-shaped knots.
#'
#' @importFrom splines interpSpline
#' @param senna A `SENNA` object
#' @param knot_df A data frame of knots
#' @return A `SENNA` object with closed islet-shaped curve axis stored
#' @export

trtype_islet <- function(senna,
                         knot_df){

  knot_df <- rbind(knot_df, knot_df[1,])
  t <- seq_along(knot_df$X1)

  if(length(t) <= 3){
    stop("To generate closed curve, at least three knots is needed.")
  }

  else{
    res <- ctype_spline(senna = senna,
                        knot_df = knot_df,
                        t = t)
  }


  xcoef <- res@CurveAxis[["fun"]][["x.coef"]]
  ycoef <- res@CurveAxis[["fun"]][["y.coef"]]
  vec1 <- c(xcoef[1, 3], ycoef[1, 3])
  vec2 <- c(xcoef[nrow(xcoef), 3], ycoef[nrow(ycoef), 3])

  res@CurveAxis[["fun"]][["t"]] <- seq_along(knot_df$X1)
  res@CurveAxis[["fun"]][["x.coef"]] <- xcoef[2:(nrow(xcoef)-1),]
  res@CurveAxis[["fun"]][["y.coef"]] <- ycoef[2:(nrow(ycoef)-1),]
  res@CurveAxis[["exd"]] <- list(vec1 = vec1,
                                 vec2 = vec2)
  res@CurveAxis[["type"]] <- c("islet", "spline", "trimmed")


  return(res)
}



#' Compute value for straight curve axis
#'
#' Computes spatial coordinate from a curve parameter using linear curve axis (CA) coefficients.
#'
#' @param val A curve parameter (scalar)
#' @param coef_mat A coefficient matrix for linear segments
#' @param t_seq A numeric vector of break points along the curve
#' @return A numeric value representing the spatial coordinate
#' @export

lin_plug <- function(val, coef_mat, t_seq) {
  if(val < 1) {
    res <- coef_mat[1, 1] * (val) + coef_mat[1, 2]
  } else if(val >= max(t_seq)) {
    res <- coef_mat[nrow(coef_mat), 1] * (val - max(t_seq)) + coef_mat[nrow(coef_mat), 2]
  } else{
    res <- coef_mat[ceiling(val), 1] * (val - floor(val)) + coef_mat[ceiling(val), 2]
  }

  return(res)
}



#' Compute value for trimmed straight curve axis
#'
#' Computes spatial coordinate from a curve parameter using linear curve axis (CA) coefficients.
#'
#' @param val A curve parameter (scalar)
#' @param coef_mat A coefficient matrix for linear segments
#' @param t_seq A numeric vector of break points along the curve
#' @return A numeric value representing the spatial coordinate
#' @export

trlin_plug <- function(val, coef_mat, t_seq) {
  res <- coef_mat[floor(val), 1] * (val - floor(val)) + coef_mat[floor(val), 2]
  return(res)
}
