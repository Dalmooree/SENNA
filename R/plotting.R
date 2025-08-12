#' @include create_senna.R
#'
NULL

#' Visualize spatial tissue layout from a SENNA object.
#'
#' Plot spatial coordinates from a `SENNA` object to visualize the tissue field-of-view (FOV).  
#' If annotation exists in `Gene$Reference$Annotation`, color spots by cell type.
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 scale_color_manual
#' @importFrom dplyr mutate
#'
#' @param senna A `SENNA` object
#' @param dot_size Numeric value for dot size (default = 1)
#' @param dot_alpha Numeric value for dot transparency (default = 1)
#' @param colors Optional named vector of colors for manual cell type coloring
#'
#' @return A `ggplot` object representing the tissue spatial layout
#' @export

TissueFOV <- function(senna,
                      dot_size = 1,
                      dot_alpha = 1,
                      colors = NULL) {

  coord <- senna@Coord[["Spatial"]]

  if(!is.null(senna@Gene[["Reference"]][["Annotation"]])){
    coord <- merge(x = coord,
                   y = data.frame(senna@Gene[["Reference"]][["Annotation"]]),
                   by = "row.names",
                   all = TRUE)
    rownames(coord) <- coord[, 1]
    coord <- coord[, 2:4]
    colnames(coord) <- c("X1", "X2", "CLS")
  }

  else{
    coord <- dplyr::mutate(coord,
                           CLS = "grey")
  }

  if(!is.null(colors)){
    coord %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(aes(X1, X2, col = CLS),
                          size = dot_size, alpha = dot_alpha) +
      ggplot2::theme_test() +
      ggplot2::scale_color_manual(values = colors)
  }

  else{
    coord %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(aes(X1, X2, col = CLS),
                          size = dot_size, alpha = dot_alpha) +
      ggplot2::theme_test()
  }
}



#' Visualize the curve axis from a SENNA object.
#'
#' Defines the generic `ShowCurve()` function for plotting the fitted curve axis  
#' along with background spatial coordinates. See method documentation for specific plotting behavior.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A SENNA object
#' @param ... Additional arguments passed to the method
#' @export
#' @seealso \code{\link[=ShowCurve,SENNA-method]{ShowCurve for SENNA}}
#'
setGeneric("ShowCurve", function(senna, ...) {
  standardGeneric("ShowCurve")
})




#' Plot the curve axis and knots over the tissue.
#'
#' Visualize the curve axis and its associated knot points overlaid on the tissue layout.  
#' If a cell-type annotation is provided via `color_reference`, the tissue background is colored accordingly.
#'
#' @importFrom methods setMethod
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna A `SENNA` object with curve axis information
#' @param order_label If `TRUE`, show index numbers next to each knot along the curve
#' @param color_reference A user-defined reference name pointing to a vector added via `AddReference()`, stored in `Gene$Reference`
#' @param bg_dot_size Size of the background tissue points
#' @param bg_dot_alpha Transparency of the background tissue points
#' @param colors Optional color specification (single color or named vector for manual mapping)
#' @param line_size Size of the curve trace line
#' @param line_color Color of the curve trace line
#' @param knots_size Size of knot points along the curve
#' @param knots_color Color of knot points
#' @param text_size Text size used for knot index labels (if `order_label = TRUE`)
#'
#' @return A `ggplot` object with curve axis and tissue context
#' @exportMethod ShowCurve
#' @seealso \code{\link[=ShowCurve]{ShowCurve}}

setMethod("ShowCurve",
          signature = list(senna = "SENNA"),
          function(senna,
                   order_label = TRUE,
                   color_reference = NULL,
                   bg_dot_size = 1,
                   bg_dot_alpha = 1,
                   colors = NULL,
                   line_size = 1,
                   line_color = 'skyblue',
                   knots_size = 1.5,
                   knots_color = 'blue',
                   text_size = 3
          ){

            griddat <- crvtrjry(senna)

            if(!is.null(color_reference)) {
              val <- senna@Coord[["Spatial"]][,c("X1", "X2")]
              val <- merge(x = val,
                           y = data.frame(senna@Gene[["Reference"]][[color_reference]]),
                           by = "row.names",
                           all = TRUE)
              rownames(val) <- val[, 1]
              val <- val[, 2:4]
              colnames(val) <- c("X1", "X2", "CLS")

              if(!is.null(colors)) {
                p <- ggplot2::ggplot() +
                  ggplot2::geom_point(aes(X1, X2, color = CLS),
                                      data = val,
                                      alpha = bg_dot_alpha,
                                      size = bg_dot_size) +
                  ggplot2::scale_color_manual(values = colors)
              }

              else{
                p <- ggplot2::ggplot() +
                  ggplot2::geom_point(aes(X1, X2, color = CLS),
                                      data = val,
                                      alpha = bg_dot_alpha,
                                      size = bg_dot_size)
              }

            }

            else{
              if(is.null(colors)){
                p <- ggplot2::ggplot(data = NULL, aes(X1, X2)) +
                  ggplot2::geom_point(data = senna@Coord[["Spatial"]],
                                      alpha = bg_dot_alpha,
                                      color = "#cccccc",
                                      size = bg_dot_size)
              } else{
                p <- ggplot2::ggplot(data = NULL, aes(X1, X2)) +
                  ggplot2::geom_point(data = senna@Coord[["Spatial"]],
                                      alpha = bg_dot_alpha,
                                      color = colors,
                                      size = bg_dot_size) 
              }
            }

            p <- p +
              ggplot2::geom_path(aes(X1, X2),
                                  data = griddat, color = line_color, size = line_size) +
              ggplot2::geom_point(aes(X1, X2),
                                  data = senna@CurveAxis[["knots"]], color = knots_color, size = knots_size) +
              ggplot2::theme_test()

            if(order_label == TRUE){
              p +
                ggrepel::geom_text_repel(aes(X1, X2),
                                         label = 1:nrow(senna@CurveAxis[["knots"]]),
                                         data = senna@CurveAxis[["knots"]],
                                         max.overlaps = nrow(senna@CurveAxis[["knots"]]),
                                         size = text_size
                )
            } else{
              p
            }
          }
)



#' Generate curve trajectory coordinates
#'
#' Generate intermediate points along the curve axis based on its type.  
#' Automatically dispatches to one of the internal trajectory generators:  
#' linear (`crvtrjrylin()`), spline (`crvtrjryspl()`), or 3-knot spline (`crvtrjry3k()`).
#'
#' @param senna A `SENNA` object containing the curve definition
#'
#' @return A data frame of 2D spatial coordinates sampled along the curve
#' @export

crvtrjry <- function(senna){

  if("straight" %in% senna@CurveAxis[["type"]]) {
    griddat <- crvtrjrylin(senna)
  }

  else if("spline" %in% senna@CurveAxis[["type"]]) {
    if (senna@CurveAxis[["ord"]] == 3){
      griddat <- crvtrjryspl(senna)
    }

    else{
      griddat <- crvtrjry3k(senna)
    }
  }

  return(griddat)
}



#' Generate spline-based curve trajectory
#'
#' Generate a dense set of spatial coordinates along a cubic spline curve  
#' interpolated through the user-defined knots. If the curve is not trimmed,  
#' linear extrapolation is added at both ends using precomputed curve parameters.
#'
#' @importFrom splines interpSpline
#'
#' @param senna A `SENNA` object containing spline curve axis information
#'
#' @return A data frame of 2D coordinates sampled along the spline trajectory
#' @export

crvtrjryspl <- function(senna){
  # Generate parameter t.
  knot.df <- senna@CurveAxis[["knots"]]
  t <- seq_along(knot.df$X1)
  Traj <- seq(min(t), max(t), length.out = (max(t) - 1) * 250)

  spx <- splines::interpSpline(knot.df$X1 ~ t)
  spy <- splines::interpSpline(knot.df$X2 ~ t)

  grid2 <- data.frame(t = Traj,
                      X1 = stats::predict(spx, Traj)$y,
                      X2 = stats::predict(spy, Traj)$y)

  if("trimmed" %in% senna@CurveAxis[["type"]]) {
    griddat <- grid2
  }
  else{
    fun <- senna@CurveAxis[["fun"]]
    xx <- fun$x.coef
    yy <- fun$y.coef
    t2 <- fun$t

    ts <- seq(min(t2), 1, length.out = 101)[1:100]
    grid1 <- data.frame(t = ts,
                        X1 = xx[1, 3] * ts + xx[1, 4],
                        X2 = yy[1, 3] * ts + yy[1, 4])

    tl <- seq(max(t), max(t2), length.out = 101)[2:101]
    grid3 <- data.frame(t = tl,
                        X1 = xx[nrow(xx), 3] * (tl - max(t)) + xx[nrow(xx), 4],
                        X2 = yy[nrow(yy), 3] * (tl - max(t)) + yy[nrow(yy), 4])


    griddat <- rbind(grid1, grid2, grid3)
  }

  return(griddat)
}



#' Generate 3-knot curve trajectory using polynomial coefficients
#'
#' Generate a piecewise curve trajectory when the curve axis is defined  
#' with exactly three knots. Polynomial interpolation is applied to the  
#' precomputed coefficients to sample dense 2D spatial coordinates.  
#' Linear extrapolation is supported if the curve is not trimmed.
#'
#' @importFrom pracma polyval
#'
#' @param senna A `SENNA` object containing a curve axis with three knots
#'
#' @return A data frame of 2D coordinates sampled along the multi-segment trajectory
#' @export

crvtrjry3k <- function(senna){

  knot.df <- senna@CurveAxis[["knots"]]
  t <- senna@CurveAxis[["fun"]][["t"]]
  tgrid <- seq(0, 1, length.out = 250)
  xcoef <- senna@CurveAxis[["fun"]][["x.coef"]]
  ycoef <- senna@CurveAxis[["fun"]][["y.coef"]]

  if("trimmed" %in% senna@CurveAxis[["type"]]) {
    tgrid <- seq(0, 1, length.out = 250)

    xval <- c(sapply(1:2, function(i) pracma::polyval(xcoef[i,], tgrid)))
    yval <- c(sapply(1:2, function(i) pracma::polyval(ycoef[i,], tgrid)))
    griddat <- data.frame(t = c(tgrid + 1, tgrid + 2),
                          X1 = xval,
                          X2 = yval)
  }
  else{
    tgrdm <- seq(min(t), 1, length.out = 250)
    tgrid <- seq(0, 1, length.out = 250)
    tgrdx <- seq(0, max(t) - 3, length.out = 250)

    xval <- c(pracma::polyval(xcoef[1,], tgrdm),
              c(sapply(2:3, function(i) pracma::polyval(xcoef[i,], tgrid))),
              pracma::polyval(xcoef[4,], tgrdx))

    yval <- c(pracma::polyval(ycoef[1,], tgrdm),
              c(sapply(2:3, function(i) pracma::polyval(ycoef[i,], tgrid))),
              pracma::polyval(ycoef[4,], tgrdx))

    griddat <- data.frame(t = c(tgrdm, tgrid + 1, tgrid + 2, tgrdx + 3),
                          X1 = xval,
                          X2 = yval)
  }

  return(griddat)
}



#' Generate straight-line curve trajectory using piecewise linear interpolation
#'
#' Compute intermediate spatial coordinates along a curve axis when  
#' it is defined as a straight path. This method uses piecewise linear  
#' interpolation of the precomputed coefficients.
#'
#' @param senna A `SENNA` object with a straight-type curve axis
#'
#' @return A data frame of 2D coordinates sampled along the straight-line trajectory
#' @export

crvtrjrylin <- function(senna){
  fun <- senna@CurveAxis[["fun"]]
  lnx <- fun$x.coef
  lny <- fun$y.coef
  t_seq <- seq_along(senna@CurveAxis[["knots"]]$X1)
  t <- fun$t
  param <- seq(min(t), max(t), length.out = (max(t_seq) + 1) * 250)
  
  if("trimmed" %in% senna@CurveAxis[["type"]]){
    griddat <- data.frame(t = param,
                          X1 = sapply(param, function(i) trlin_plug(i, lnx, t_seq)),
                          X2 = sapply(param, function(i) trlin_plug(i, lny, t_seq)))
  } else{
    griddat <- data.frame(t = param,
                          X1 = sapply(param, function(i) lin_plug(i, lnx, t_seq)),
                          X2 = sapply(param, function(i) lin_plug(i, lny, t_seq)))
  }

  return(griddat)
}



#' Plot the generated curve axis
#'
#' Visualize the fitted curve axis over the spatial layout.  
#' Depending on whether the curve is trimmed or not, this function delegates  
#' to either `covalwtr()` (trimmed) or `covalwotr()` (untrimmed).
#'
#' @param senna A `SENNA` object containing a fitted curve axis
#' @param seed Random seed for plotting consistency
#' @param bg_dot_size Size of background tissue points
#' @param bg_dot_alpha Transparency of background tissue points
#' @param lead_size Line width of curve guidance overlay
#' @param lead_type Line type for curve guidance (e.g., `"longdash"`, `"solid"`)
#' @param line_size Size of the curve axis trace
#'
#' @return A `ggplot` object showing the curve axis on tissue coordinates
#' @export

CurveParamValid <- function(senna,
                            seed = 1,
                            bg_dot_size = 0.7,
                            bg_dot_alpha = 0.2,
                            lead_size = 1,
                            lead_type = "longdash",
                            line_size = 1) {
  if(!"trimmed" %in% senna@CurveAxis[["type"]]) {
    covalwotr(senna, seed, bg_dot_size, bg_dot_alpha, lead_size, lead_type, line_size)
  }
  else if("trimmed" %in% senna@CurveAxis[["type"]]){
    covalwtr(senna, seed, bg_dot_size, bg_dot_alpha, lead_size, lead_type, line_size)
  }
  else{
    stop("Covariate warning : Invalid curve")
  }
}



#' Plot full curve with random lead lines
#'
#' Visualize the fitted curve axis overlaid on tissue coordinates.  
#' Randomly selects several spots and connects them to the corresponding curve projections  
#' to illustrate alignment between spatial coordinates and the curve.
#'
#' @importFrom magrittr %>%
#' @importFrom BiocGenerics rbind
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 theme
#'
#' @param senna A `SENNA` object with a fitted curve axis
#' @param seed Random seed for reproducible spot sampling
#' @param bg.dot.size Size of background dots representing tissue coordinates
#' @param alpha Transparency of background dots
#' @param lead Line width for lead lines connecting spot to curve
#' @param type Line type (e.g., "longdash") for lead lines
#' @param line Size of the curve axis trace
#'
#' @return A `ggplot` object visualizing the curve axis and alignment
#' @export

covalwotr <- function(senna,
                      seed,
                      bg.dot.size,
                      alpha,
                      lead,
                      type,
                      line) {
  dat <- senna@Coord[["Spatial"]]
  fun <- senna@CurveAxis[["fun"]]
  griddat <- crvtrjry(senna)

  set.seed(seed)
  ind <- sample(1:nrow(dat), 6)

  if("spline" %in% senna@CurveAxis[["type"]]) {

    d <- lapply(ind, function(i){
      subdat <- dat[i, ] %>%
        dplyr::mutate(index = i)

      subdat2 <- data.frame(X1 = plug_coef(subdat$t, fun$x.coef),
                            X2 = plug_coef(subdat$t, fun$y.coef),
                            index = i)

      subdat <- subdat %>% dplyr::select(X1, X2, index)

      sd <- BiocGenerics::rbind(subdat, subdat2)
      return(sd)
    })

  } else if("straight" %in% senna@CurveAxis[["type"]]) {

    d <- lapply(ind, function(i){
      subdat <- dat[i, ] %>%
        dplyr::mutate(index = i)

      subdat2 <- data.frame(X1 = plug_coef_lin(subdat$t, fun$x.coef),
                            X2 = plug_coef_lin(subdat$t, fun$y.coef),
                            index = i)

      subdat <- subdat %>% dplyr::select(X1, X2, index)

      sd = BiocGenerics::rbind(subdat, subdat2)
      return(sd)
    })

  }


  ggplot2::ggplot()+
    ggplot2::geom_point(aes(X1, X2), data = dat, col = 'grey', alpha = alpha, size = bg.dot.size) +
    ggplot2::geom_path(aes(X1, X2), data = griddat, color = 'skyblue', linewidth = line) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#b14743", size = lead, data = d[[1]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#44729d", size = lead, data = d[[2]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#d48640", size = lead, data = d[[3]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#8e73ae", size = lead, data = d[[4]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#539045", size = lead, data = d[[5]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#7e5d55", size = lead, data = d[[6]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#b14743", fill = 'grey',
                        data = d[[1]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#44729d", fill = 'grey',
                        data = d[[2]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#d48640", fill = 'grey',
                        data = d[[3]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#8e73ae", fill = 'grey',
                        data = d[[4]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#539045", fill = 'grey',
                        data = d[[5]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#7e5d55", fill = 'grey',
                        data = d[[6]]) +
    ggplot2::theme_test() +
    ggplot2::theme(legend.position = "none")
}



#' Plot trimmed curve with lead lines
#'
#' Visualize a trimmed curve axis and illustrate alignment between selected tissue spots  
#' and their projected positions along the curve. Spots with missing curve parameter (`t`)  
#' are shown in a lighter color to distinguish them from valid points.
#'
#' @importFrom magrittr %>%
#' @importFrom BiocGenerics rbind
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 theme
#'
#' @param senna A `SENNA` object containing a trimmed curve
#' @param seed Random seed for reproducible spot sampling
#' @param bg.dot.size Size of background dots representing tissue coordinates
#' @param alpha Transparency level for background and foreground dots
#' @param lead Line width for the lead lines connecting tissue spots to the curve
#' @param type Line type (e.g., "longdash") for the lead lines
#' @param line Size of the curve axis trace
#'
#' @return A `ggplot` object visualizing trimmed curve behavior and spot alignment
#' @export

covalwtr <- function(senna,
                     seed,
                     bg.dot.size,
                     alpha,
                     lead,
                     type,
                     line) {
  dat <- subset(senna@Coord[["Spatial"]], !is.na(t))
  dat1 <- subset(senna@Coord[["Spatial"]], is.na(t))
  fun <- senna@CurveAxis[["fun"]]
  griddat <- crvtrjry(senna)

  set.seed(seed)
  ind <- sample(1:nrow(dat), 6)

  if("spline" %in% senna@CurveAxis[["type"]]) {

    d <- lapply(ind, function(i){
      subdat <- dat[i, ] %>%
        dplyr::mutate(index = i)

      subdat2 <- data.frame(X1 = trplug_coef(subdat$t, fun$x.coef),
                            X2 = trplug_coef(subdat$t, fun$y.coef),
                            index = i)

      subdat <- subdat %>% dplyr::select(X1, X2, index)

      sd = BiocGenerics::rbind(subdat, subdat2)
      return(sd)
    })

  } else if("straight" %in% senna@CurveAxis[["type"]]) {

    d <- lapply(ind, function(i){
      subdat <- dat[i, ] %>%
        dplyr::mutate(index = i)

      subdat2 <- data.frame(X1 = trplug_coef_lin(subdat$t, fun$x.coef),
                            X2 = trplug_coef_lin(subdat$t, fun$y.coef),
                            index = i)

      subdat <- subdat %>% dplyr::select(X1, X2, index)

      sd = BiocGenerics::rbind(subdat, subdat2)
      return(sd)
    })

  }


  ggplot2::ggplot()+
    ggplot2::geom_point(aes(X1, X2), data = dat1, col = 'lightgrey', alpha = alpha, size = bg.dot.size) +
    ggplot2::geom_point(aes(X1, X2), data = dat, col = 'darkgrey', alpha = alpha, size = bg.dot.size) +
    ggplot2::geom_path(aes(X1, X2), data = griddat, color = 'skyblue', linewidth = line) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#b14743", size = lead, data = d[[1]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#44729d", size = lead, data = d[[2]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#d48640", size = lead, data = d[[3]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#8e73ae", size = lead, data = d[[4]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#539045", size = lead, data = d[[5]]) +
    ggplot2::geom_line(aes(X1, X2), linetype = type,
                       col = "#7e5d55", size = lead, data = d[[6]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#b14743", fill = 'grey',
                        data = d[[1]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#44729d", fill = 'grey',
                        data = d[[2]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#d48640", fill = 'grey',
                        data = d[[3]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#8e73ae", fill = 'grey',
                        data = d[[4]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#539045", fill = 'grey',
                        data = d[[5]]) +
    ggplot2::geom_point(aes(X1, X2),
                        shape = 21, size = bg.dot.size * 1.2,
                        col = "#7e5d55", fill = 'grey',
                        data = d[[6]]) +
    ggplot2::theme_test() +
    ggplot2::theme(legend.position = "none")
}



#' Volcano plot for progression SVGs
#'
#' Define the generic function `ProgVolPlot()` to visualize progression SVGs  
#' from either a `SENNA` or `mSENNA` object. See method implementations for details.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments passed to the corresponding method
#'
#' @return A `ggplot` object representing a volcano plot of progression SVGs
#' @export
#' @seealso \code{\link[=ProgVolPlot,SENNA-method]{ProgVolPlot for SENNA}},  
#'   \code{\link[=ProgVolPlot,msR-method]{ProgVolPlot for msR}}


setGeneric("ProgVolPlot", function(senna, ...) {
  standardGeneric("ProgVolPlot")
})



#' Create volcano plot for progression SVGs (SENNA method)
#'
#' Generates a volcano plot summarizing progression SVGs identified in a `SENNA` object.  
#' Significant genes are separated into positive and negative gradients based on sign and magnitude.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr arrange
#' @importFrom dplyr setdiff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggrepel geom_label_repel
#'
#' @param senna A `SENNA` object containing progression SVGs
#' @param FDR_level FDR threshold for significance cutoff
#' @param grad_cutoff Threshold for absolute value of gradient coefficients
#' @param dot_size Size of points in the plot
#' @param dot_alpha Transparency of points
#' @param positive Color for positively associated genes
#' @param negative Color for negatively associated genes
#' @param cutoff_col Color for cutoff lines (FDR and gradient)
#' @param nrepel Number of gene labels to display using text repelling
#' @param text_size Size of gene labels shown by `nrepel`
#' @param yfix Logical, if `TRUE` fixes upper y-axis limit to improve comparability
#'
#' @return A `ggplot` object representing the volcano plot of progression SVGs
#' @exportMethod ProgVolPlot

setMethod("ProgVolPlot",
          signature(senna = "SENNA"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0,
                   dot_size = 0.5,
                   dot_alpha = 0.5,
                   positive = "#ec1b24",
                   negative = "#1b74bc",
                   cutoff_col = "#888888",
                   nrepel = 0L,
                   text_size = 3,
                   yfix = FALSE){

            report <- senna@Gene[["P.SVGs"]][["Report"]]
            possig <- subset(report, Gene %in% senna@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]])
            negsig <- subset(report, Gene %in% senna@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]])
            nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

            if(yfix){
              p <- ggplot2::ggplot() +
                ggplot2::ylim(c(NA, 301)) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                     col = positive, alpha = dot_alpha,
                                     size = dot_size, data = possig) +
                ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                     col = negative, alpha = dot_alpha,
                                     size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_test() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            } else{
              p <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                     col = positive, alpha = dot_alpha,
                                     size = dot_size, data = possig) +
                ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                     col = negative, alpha = dot_alpha,
                                     size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_test() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            }

            if(grad_cutoff != 0){
              p <- p +
                ggplot2::geom_vline(xintercept = sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col) +
                ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col)
            }

            if(nrepel > 0){
              ntop <- dplyr::arrange(report, log(adj.p))[1:nrepel, ]
              p <- p +
                ggrepel::geom_text_repel(
                  aes(Coefficients / sqrt(abs(Coefficients)),
                      -log10(adj.p + 1e-300),
                      label = Gene),
                  size = text_size,
                  data = ntop)
            }

            p
          })



#' Volcano plot for progression SVGs (msR method)
#'
#' Generates a volcano plot for progression-related SVGs stored in an `msR` object.  
#' Positive and negative gradient genes are plotted with distinct colors, and optional labels can highlight selected genes.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr arrange
#' @importFrom dplyr setdiff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna An `msR` object containing progression SVG analysis results
#' @param FDR_level FDR threshold for significance cutoff
#' @param grad_cutoff Threshold for the absolute value of gradient coefficients
#' @param dot_size Size of the points in the plot
#' @param dot_alpha Transparency of points
#' @param positive Color for positively associated genes
#' @param negative Color for negatively associated genes
#' @param cutoff_col Color used for cutoff lines (FDR or gradient)
#' @param nrepel Number of top genes to label using text repelling (default = 0)
#' @param text_size Text size for repelled gene labels
#' @param yfix Logical. If `TRUE`, fixes the upper y-axis limit for comparability
#'
#' @return A `ggplot` object representing the volcano plot of progression SVGs
#' @exportMethod ProgVolPlot

setMethod("ProgVolPlot",
          signature(senna = "msR"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0,
                   dot_size = 0.5,
                   dot_alpha = 0.5,
                   positive = "#ec1b24",
                   negative = "#1b74bc",
                   cutoff_col = "#888888",
                   nrepel = 0L,
                   text_size = 3,
                   yfix = FALSE){

            report <- senna@Report
            possig <- subset(report, Gene %in% senna@Variable_gene[["positive"]])
            negsig <- subset(report, Gene %in% senna@Variable_gene[["negative"]])
            nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

            if(yfix){
              p <- ggplot2::ggplot() +
                ggplot2::ylim(c(NA, 301)) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = positive, alpha = dot_alpha,
                                    size = dot_size, data = possig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = negative, alpha = dot_alpha,
                                    size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_light() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            } else{
              p <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = positive, alpha = dot_alpha,
                                    size = dot_size, data = possig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = negative, alpha = dot_alpha,
                                    size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_light() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            }

            if(grad_cutoff != 0){
              p <- p +
                ggplot2::geom_vline(xintercept = sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col) +
                ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col)
            }

            if(nrepel > 0){
              ntop <- dplyr::arrange(report, log(adj.p))[1:nrepel, ]
              p <- p +
                ggrepel::geom_text_repel(
                  aes(Coefficients / sqrt(abs(Coefficients)),
                      -log10(adj.p + 1e-300),
                      label = Gene),
                  size = text_size,
                  data = ntop)
            }

            p
          })



#' Plot interval for progression analysis
#'
#' Visualize spatial spots within a specified interval from the curve axis.  
#' Colors represent either binary box weights or continuous Gaussian weights  
#' depending on the `gaussian` option.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 theme
#' @importFrom stats dnorm
#'
#' @param senna A `SENNA` object with spatial coordinates and curve axis
#' @param interval A numeric distance threshold used for spot highlighting
#' @param gaussian If `TRUE`, show continuous Gaussian weights; otherwise, show binary box weights
#' @param dot_size Size of tissue spots in the plot
#' @param dot_alpha Transparency of tissue spots
#'
#' @return A `ggplot` object showing the interval-weighted tissue layout
#' @export

ScanInterval = function(senna,
                        interval = 1,
                        gaussian = FALSE,
                        dot_size = 1,
                        dot_alpha = 1) {

  if(gaussian == FALSE){
    coord <- senna@Coord[["Spatial"]] %>%
      dplyr::mutate(box = ifelse(distance <= interval, 1, 0))

    coord %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(aes(X1, X2, col = box), size = dot_size, alpha = dot_alpha) +
      ggplot2::scale_colour_gradient(high = "skyblue", low = "darkgrey", na.value = "lightgrey") +
      ggplot2::theme_test() +
      ggplot2::theme(legend.position = "none")

  } else if(gaussian == TRUE){
    coord <- senna@Coord[["Spatial"]] %>%
      dplyr::mutate(gwei = stats::dnorm(distance) / stats::dnorm(0))

    ggplot2::ggplot() +
      ggplot2::geom_point(aes(X1, X2), col = "lightgrey", data = coord, alpha = dot_alpha, size = dot_size) +
      ggplot2::geom_point(aes(X1, X2, col = gwei), data = dplyr::filter(coord, distance <= interval), alpha = dot_alpha, size = dot_size) +
      ggplot2::scale_colour_gradient(low ='lightgrey', high = "darkblue") +
      ggplot2::theme_test() +
      ggplot2::theme(legend.position = "none")
  }
}



#' Volcano plot for regionation SVGs
#'
#' Define generic function `RegionVolPlot`.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments passed to the method
#'
#' @return A `ggplot` object representing the volcano plot
#' @export
#' @seealso \code{\link[=RegionVolPlot,SENNA-method]{RegionVolPlot for SENNA}},
#'   \code{\link[=RegionVolPlot,msR-method]{RegionVolPlot for msR}}

setGeneric("RegionVolPlot", function(senna, ...) {
  standardGeneric("RegionVolPlot")
})



#' Volcano plot for regionation SVGs
#'
#' Visualize spatially variable genes (SVGs) identified from regionation analysis  
#' using a volcano plot. Shows effect size versus statistical significance.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr arrange
#' @importFrom dplyr setdiff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna A `SENNA` object
#' @param FDR_level False discovery rate threshold for significance
#' @param grad_cutoff Threshold for absolute value of gradients (effect size)
#' @param dot_size Size of dots in the volcano plot
#' @param dot_alpha Transparency of dots
#' @param positive Color for positive SVGs (upregulated)
#' @param negative Color for negative SVGs (downregulated)
#' @param cutoff_col Color for cutoff lines
#' @param nrepel Number of top genes to label
#' @param text_size Size of gene labels
#' @param yfix Whether to fix the vertical axis limit (default = FALSE)
#'
#' @return A `ggplot` object representing the volcano plot
#' @exportMethod RegionVolPlot

#'
setMethod("RegionVolPlot",
          signature(senna = "SENNA"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0,
                   dot_size = 0.5,
                   dot_alpha = 0.5,
                   positive = "#98499a",
                   negative = "#67b665",
                   cutoff_col = "#888888",
                   nrepel = 0L,
                   text_size = 3,
                   yfix = FALSE
                   ){

            report <- senna@Gene[["R.SVGs"]][["Report"]]
            possig <- subset(report, Gene %in% senna@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]])
            negsig <- subset(report, Gene %in% senna@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])
            nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

            if(yfix){
              p <- ggplot2::ggplot() +
                ggplot2::ylim(c(NA, 301)) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = positive, alpha = dot_alpha,
                                    size = dot_size, data = possig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = negative, alpha = dot_alpha,
                                    size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_light() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
              } else{
                p <- ggplot2::ggplot() +
                  ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                      col = 'gray', alpha = dot_alpha,
                                      size = dot_size, data = nonsig) +
                  ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                      col = positive, alpha = dot_alpha,
                                      size = dot_size, data = possig) +
                  ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                      col = negative, alpha = dot_alpha,
                                      size = dot_size, data = negsig) +
                  ggplot2::labs(x = "Gradients (Scaled)",
                                y = "Adjusted p-value, -log") +
                  ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                  ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                  ggplot2::theme_light() +
                  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                 axis.ticks.y = ggplot2::element_blank())
                }

            if(grad_cutoff != 0){
              p <- p +
                ggplot2::geom_vline(xintercept = sqrt(grad_cutoff),
                                    linetype = "dotted", alpha = 0.7,
                                    linewidth = 1.3, col = cutoff_col) +
                ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff),
                                    linetype = "dotted", alpha = 0.7,
                                    linewidth = 1.3, col = cutoff_col)
            }

            if(nrepel > 0){
              ntop <- dplyr::arrange(report, log(adj.p))[1:nrepel, ]
              p <- p +
                ggrepel::geom_text_repel(
                  aes(Coefficients / sqrt(abs(Coefficients)),
                      -log10(adj.p + 1e-300),
                      label = Gene),
                  size = text_size,
                  data = ntop)
            }

            p
            })



#' Volcano plot for Regionation SVGs
#'
#' Visualize spatially variable genes (SVGs) identified from regionation analysis  
#' using a volcano plot. Shows effect size versus statistical significance.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr arrange
#' @importFrom dplyr setdiff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna An `msR` object (output from regionation peak detection)
#' @param FDR_level False discovery rate threshold for significance
#' @param grad_cutoff Threshold for absolute value of gradients (effect size)
#' @param dot_size Size of dots in the volcano plot
#' @param dot_alpha Transparency of dots
#' @param positive Color for positive SVGs (upregulated)
#' @param negative Color for negative SVGs (downregulated)
#' @param cutoff_col Color for significance threshold lines
#' @param nrepel Number of top genes to label using `geom_text_repel`
#' @param text_size Size of the text for gene labels
#' @param yfix If `TRUE`, fix the vertical axis limit
#'
#' @return A `ggplot` object representing the volcano plot
#' @exportMethod RegionVolPlot

setMethod("RegionVolPlot",
          signature(senna = "msR"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0,
                   dot_size = 0.5,
                   dot_alpha = 0.5,
                   positive = "#98499a",
                   negative = "#67b665",
                   cutoff_col = "#888888",
                   nrepel = 0L,
                   text_size = 3,
                   yfix = FALSE){

            report <- senna@Report
            possig <- subset(report, Gene %in% senna@Variable_gene[["positive"]])
            negsig <- subset(report, Gene %in% senna@Variable_gene[["negative"]])
            nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

            if(yfix){
              p <- ggplot2::ggplot() +
                ggplot2::ylim(c(NA, 301)) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = positive, alpha = dot_alpha,
                                    size = dot_size, data = possig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = negative, alpha = dot_alpha,
                                    size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_light() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            } else{
              p <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = 'gray', alpha = dot_alpha,
                                    size = dot_size, data = nonsig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = positive, alpha = dot_alpha,
                                    size = dot_size, data = possig) +
                ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                    col = negative, alpha = dot_alpha,
                                    size = dot_size, data = negsig) +
                ggplot2::labs(x = "Gradients (Scaled)",
                              y = "Adjusted p-value, -log") +
                ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                ggplot2::theme_light() +
                ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
            }

            if(grad_cutoff != 0){
              p <- p +
                ggplot2::geom_vline(xintercept = sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col) +
                ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col)
            }

            if(nrepel > 0){
              ntop <- arrange(report, log(adj.p))[1:nrepel, ]
              p <- p +
                ggrepel::geom_text_repel(
                  aes(Coefficients / sqrt(abs(Coefficients)),
                      -log10(adj.p + 1e-300),
                      label = Gene),
                  size = text_size,
                  data = ntop)
            }

            p
          })



#' generics for CSD plot
#'
#' Define generic function `ShowCSDistance()` for plotting C-S distance from the curve axis.  
#' Methods are implemented for both `SENNA` and `mSENNA` objects.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments passed to the method
#'
#' @return A `ggplot` object visualizing C-S distance
#' @export
#' @seealso \code{\link[=ShowCSDistance,SENNA-method]{ShowCSDistance for SENNA}},  
#'   \code{\link[=ShowCSDistance,mSENNA-method]{ShowCSDistance for mSENNA}}

setGeneric("ShowCSDistance", function(senna,
                                      ...) {
  standardGeneric("ShowCSDistance")
})



#' Plot C-S distance along the curve axis
#'
#' Visualize each spot’s signed or absolute distance from the curve axis.  
#' Optionally highlights one side of the region by showing the other in background color.
#'
#' @importFrom methods setMethod
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom ggplot2 theme_test
#'
#' @param senna A `SENNA` object
#' @param dot_size Size of spot dots
#' @param dot_alpha Transparency of spot dots
#' @param line_size Size of the curve axis trace
#' @param high Color for high (positive) distance
#' @param medium Color for zero or near-zero distance
#' @param low Color for low (negative) distance
#' @param bgcol Color for background spots (if `direction` is specified)
#' @param direction Optional region ID to highlight (`1` or `-1`).  
#'   Spots not in this region are shown in background color.
#' @param set_midpoint The midpoint value used for `scale_colour_gradient2()` (only when `direction = NULL`)
#'
#' @return A `ggplot` object showing spatial distribution of C-S distances
#' @export

setMethod("ShowCSDistance",
          signature = list(senna = "SENNA"),
          function(senna,
                   dot_size = 1,
                   dot_alpha = 1,
                   line_size = 1,
                   high = "#98499a",
                   medium = "#bbbbbb",
                   low = "#67b665",
                   bgcol = "#dddddd",
                   direction = NULL,
                   set_midpoint = 0) {
            dat <- senna@Coord[["Spatial"]]
            griddat <- crvtrjry(senna)
            
            if(!is.null(direction)){
              valsub <- split(dat, dat[["region"]] == direction)
              valsub1 <- valsub[["TRUE"]]
              valsub2 <- valsub[["FALSE"]]
              
              g <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(X1, X2, col = abs(distance)), 
                                    alpha = dot_alpha, 
                                    size = dot_size,
                                    data = valsub1) +
                ggplot2::scale_colour_gradient(low = low, 
                                               high = high) +
                ggplot2::geom_point(aes(X1, X2),
                                    alpha = dot_alpha,
                                    size = dot_size, 
                                    col = bgcol,
                                    data = valsub2) +
                ggplot2::geom_path(aes(X1, X2), 
                                   data = griddat, 
                                   color = "black",
                                   size = line_size) +
                ggplot2::theme_test()
            } else{
              g <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(X1, X2, col = distance), 
                                    data = dat, 
                                    alpha = dot_alpha, 
                                    size = dot_size) +
                ggplot2::scale_colour_gradient2(high = high, 
                                                mid = medium, 
                                                low = low, 
                                                midpoint = set_midpoint) +
                ggplot2::geom_path(aes(X1, X2), 
                                   data = griddat, 
                                   color = "black", 
                                   size = line_size)
            }
            return(g)
})



#' Plot C-S distance for mSENNA object
#'
#' Visualize absolute C-S distances from multiple islet curve axes.  
#' Highlights one region (`direction`) while graying out the rest.  
#' Curve traces for each embedded SENNA object are overlaid on the background.
#'
#' @importFrom methods setMethod
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_gradient2
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom ggplot2 theme_test
#'
#' @param senna A `mSENNA` object containing multiple `SENNA` sub-objects
#' @param dot_size Size of spot dots
#' @param dot_alpha Transparency of spot dots
#' @param line_size Size of the curve trace overlay
#' @param high Color for high absolute distances
#' @param low Color for low absolute distances
#' @param bgcol Color for background (non-highlighted) region
#' @param direction Region direction to highlight (`1` or `-1`)
#' @param set_midpoint Ignored for `mSENNA` method (included for compatibility)
#'
#' @return A `ggplot` object showing absolute C-S distances with optional region highlighting
#' @export

setMethod("ShowCSDistance",
          signature = list(senna = "mSENNA"),
          function(senna,
                   dot_size = 1,
                   dot_alpha = 1,
                   line_size = 0.3,
                   high = "#ca6804",
                   low = "#dec9b3",
                   bgcol = "#dddddd",
                   direction = NULL,
                   set_midpoint = 0) {
            dat <- senna@SENNA[[1]]@Coord[["Spatial"]]
            
            if(!is.null(direction)){
              rv <- lapply(senna@SENNA, 
                           function(x) dplyr::select(x@Coord[["Spatial"]], 
                                                     "region"))
              rv <- do.call(cbind, rv)
              rv$regn <- apply(rv, 1, prod)
              rv <- dplyr::select(rv, regn)
              
              dat <- cbind(dat, rv)
              valsub <- split(dat, dat[["regn"]] == direction)
              valsub1 <- valsub[["TRUE"]]
              valsub2 <- valsub[["FALSE"]]
              
              g <- ggplot2::ggplot() +
                ggplot2::geom_point(aes(X1, X2, col = abs(distance)), 
                                    alpha = dot_alpha, 
                                    size = dot_size,
                                    data = valsub1) +
                ggplot2::scale_colour_gradient(low = low, 
                                               high = high) +
                ggplot2::geom_point(aes(X1, X2),
                                    alpha = dot_alpha,
                                    size = dot_size, 
                                    col = bgcol,
                                    data = valsub2) +
                ggplot2::theme_test()
              
              for(l in 1:length(senna@SENNA)){
                g <- g + 
                  ggplot2::geom_point(aes(X1, X2), 
                                      data = crvtrjry(senna@SENNA[[l]]), 
                                      color = "black", 
                                      size = line_size)
              }
              
            } else {
              stop("For mSENNA obejct, only islet curve axes are supported.")}
            
            return(g)
          })



#' Plot region
#'
#' Visualize region labels from a `SENNA` object over spatial coordinates.  
#' Optionally colors regions with custom palette and overlays the curve axis.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#'
#' @param senna A `SENNA` object
#' @param dot_size Size of the tissue dots
#' @param dot_alpha Transparency of the tissue dots
#' @param dot_cols Optional color palette for region labels
#' @param line_size Size of the curve line
#' @param line_col Color of the curve line
#' @param legend Whether to show the legend
#'
#' @return A `ggplot` object with region assignments overlaid
#' @export

ShowRegions <- function(senna,
                        dot_size = 1,
                        dot_alpha = 1,
                        dot_cols = NULL,
                        line_size = 1,
                        line_col = "black",
                        legend = TRUE) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(aes(X1, X2, col = factor(region)),
                        data = senna@Coord[["Spatial"]],
                        alpha = dot_alpha,
                        size = dot_size) +
    ggplot2::geom_path(aes(X1, X2),
                        data = crvtrjry(senna),
                        color = line_col, linewidth = line_size)

  if(!is.null(dot_cols)){
    p <- p +
      ggplot2::scale_color_manual(values = dot_cols) +
      theme_test()
  } else{
    p <- p + theme_test()
  }

  if(legend){
    p <- p + theme(axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
      labs(color = "Region")
  } else{
    p <- p + theme(legend.position = "none",
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
      labs(color = "Region")
  }

  p
}



#' Tissue feature plot generic
#'
#' Define the generic function for `TissueFeaturePlot()`, which visualizes  
#' gene expression over spatial coordinates in `SENNA` or `mSENNA` objects.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param gene Name of gene to plot
#' @param ... Additional arguments passed to method implementations
#'
#' @return A `ggplot` object showing spatial gene expression
#' @export
#' @seealso \code{\link{TissueFeaturePlot,SENNA-method}}, \code{\link{TissueFeaturePlot,mSENNA-method}}

setGeneric("TissueFeaturePlot", function(senna,
                                         gene,
                                         ...) {
  standardGeneric("TissueFeaturePlot")
})



#' Plot gene expression on tissue
#'
#' Visualize the expression of a specific gene on the spatial coordinate system  
#' of a `SENNA` object. Supports filtering by activation status, region, or  
#' proximity to curve axis (C-S distance), and optional overlay of the curve axis.
#'
#' @importFrom methods setMethod
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 labs
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom stats complete.cases
#' @importFrom tidyselect all_of
#'
#' @param senna A `SENNA` object
#' @param gene Name of a gene to visualize
#' @param interval Optional numeric distance threshold for inclusion based on C-S distance
#' @param active Whether to restrict to activated cells (default = FALSE)
#' @param direction Optional region direction filter (e.g., 1 or -1)
#' @param curve_axes If TRUE, overlays the curve axis (default = FALSE)
#' @param low Color for low expression values
#' @param high Color for high expression values
#' @param size Size of the gene expression points
#' @param line_size Size of the curve axis overlay
#' @param alpha Transparency level for points (0–1)
#'
#' @return A `ggplot` object showing gene expression over the tissue
#' @export
#' @seealso \code{\link{TissueFeaturePlot}}

setMethod("TissueFeaturePlot",
          signature = list(senna = "SENNA"),
          function(senna,
                   gene,
                   interval = NULL,
                   active = FALSE,
                   direction = NULL,
                   curve_axes = FALSE,
                   low = "#cfcfcf",
                   high = "#800000",
                   size = 1,
                   line_size = 1,
                   alpha = 1) {
            val <- senna@Coord[["Spatial"]]
            val <- merge(val,
                         dplyr::select(senna@Gene[["Spatial"]],
                                       tidyselect::all_of(gene)),
                         by = 0)
            val <- dplyr::mutate(val, gcount = tidyselect::all_of(gene))
            val <- dplyr::mutate(val, gcount = get(gene))

            if(!is.null(interval)){
              val$gcount <- ifelse(val$distance <= interval,
                                   val$gcount, NA)
              valsub <- split(val, stats::complete.cases(val))
              valsub1 <- valsub[["TRUE"]]
              valsub1 <- dplyr::mutate(valsub1, gcount = valsub1[[gene]])
              valsub2 <- valsub[["FALSE"]]
            } else{
              if(!is.null(direction)){
                valsub <- split(val, val[["region"]] == direction)
                valsub1 <- valsub[["TRUE"]]
                valsub1 <- dplyr::mutate(valsub1,
                                         gcount = valsub1[[gene]])
                valsub2 <- valsub[["FALSE"]]
              } else{
                valsub1 <- val
                valsub2 <- NULL
              }
            }

            if(active){
              idx <- valsub1[,1] %in% senna@Coord[["Activated"]][["cells"]]
              valsub2 <- rbind(valsub2, split(valsub1, idx)[["FALSE"]])
              valsub1 <- split(valsub1, idx)[["TRUE"]]
            }

            if(!is.null(valsub2)){
              g <- ggplot2::ggplot() +
                ggplot2::geom_point(ggplot2::aes(X1, X2, colour = gcount),
                                    size = size,
                                    alpha = alpha,
                                    data = valsub1) +
                ggplot2::scale_color_gradient(
                  low = low,
                  high = high
                ) +
                ggplot2::labs(color = paste(gene)) +
                ggplot2::geom_point(ggplot2::aes(X1, X2),
                                    color = "#dfdfdf",
                                    size = size,
                                    alpha = alpha,
                                    data = valsub2) +
                ggplot2::theme_test()
            } else{
              g <- ggplot2::ggplot() +
                ggplot2::geom_point(ggplot2::aes(X1, X2, colour = gcount),
                                    size = size,
                                    alpha = alpha,
                                    data = valsub1) +
                ggplot2::scale_color_gradient(
                  low = low,
                  high = high
                ) +
                ggplot2::labs(color = paste(gene)) +
                ggplot2::theme_test()
            }
            
            if(curve_axes){
              g <- g + 
                ggplot2::geom_path(aes(X1, X2), 
                                   data = crvtrjry(senna),             
                                   color = "black", 
                                   size = line_size)}
            return(g)
          })



#' Plot gene expression on tissue
#'
#' Visualize the expression of a specific gene on the spatial layout  
#' of an `mSENNA` object. Filters spots based on region direction or activation status,  
#' and optionally overlays curve axes.
#'
#' @importFrom methods setMethod
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 labs
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom stats complete.cases
#' @importFrom tidyselect all_of
#'
#' @param senna An `mSENNA` object
#' @param gene Name of a gene to visualize
#' @param active If `TRUE`, highlight only activated cells
#' @param direction Region direction to subset spots (`1` or `-1`)
#' @param curve_axes If `TRUE`, overlay curve axes on the tissue layout
#' @param low Color for low expression values
#' @param high Color for high expression values
#' @param size Size of the gene expression points
#' @param line_size Size of curve axis overlay
#' @param alpha Transparency level for points
#'
#' @return A `ggplot` object showing gene expression over the tissue
#' @export
#' @seealso \code{\link{TissueFeaturePlot}}

setMethod("TissueFeaturePlot",
          signature = list(senna = "mSENNA"),
          function(senna,
                   gene,
                   active = FALSE,
                   direction = NULL,
                   curve_axes = FALSE,
                   low = "#cfcfcf",
                   high = "#800000",
                   size = 1,
                   line_size = 1,
                   alpha = 1) {
            val <- senna@SENNA[[1]]@Coord[["Spatial"]]
            valg <- senna@SENNA[[1]]@Gene[["Spatial"]]
            
            if(!is.null(direction)){
              rv <- lapply(senna@SENNA, 
                           function(x) dplyr::select(x@Coord[["Spatial"]], 
                                                     "region"))
              rv <- do.call(cbind, rv)
              rv$regn <- apply(rv, 1, prod)
              rv <- dplyr::select(rv, regn)
            }
            
            val<- cbind(val, rv)
            val <- merge(val,
                         dplyr::select(valg,
                                       tidyselect::all_of(gene)),
                         by = 0)
            colnames(val) <- c(colnames(val)[-ncol(val)], "gcount")
            
            if(!is.null(direction)){
              valsub <- split(val, val[["regn"]] == direction)
              valsub1 <- valsub[["TRUE"]]
              valsub2 <- valsub[["FALSE"]]
            } else{
              stop("Error: mSENNA only supports islet curve axes. Set direction to either 1 or -1.")            
            }
            
            if(active){
              idx <- valsub1[,1] %in% senna@SENNA[[1]]@Coord[["Activated"]][["cells"]]
              valsub2 <- rbind(valsub2, split(valsub1, idx)[["FALSE"]])
              valsub1 <- split(valsub1, idx)[["TRUE"]]
            }
            
            g <- ggplot2::ggplot() +
              ggplot2::geom_point(ggplot2::aes(X1, X2, colour = gcount),
                                  size = size,
                                  alpha = alpha,
                                  data = valsub1) +
              ggplot2::scale_color_gradient(
                low = low,
                high = high) +
              ggplot2::labs(color = paste(gene)) +
              ggplot2::geom_point(ggplot2::aes(X1, X2),
                                  color = "#dfdfdf",
                                  size = size,
                                  alpha = alpha,
                                  data = valsub2) +
              ggplot2::theme_test()
            
            if(curve_axes){
              for(l in 1:length(senna@SENNA)){
                g <- g + 
                  ggplot2::geom_path(aes(X1, X2), 
                                     data = crvtrjry(senna@SENNA[[l]]), 
                                     color = "black", 
                                     size = line_size)
              }
            }
            
            return(g)
            })



#' Sketch specific pattern for each SVG
#'
#' Define generic function `SketchPatterns`.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA`, `mSENNA`, or `msR` object
#' @param ... Additional arguments passed to methods
#'
#' @export
#' @seealso \code{\link[=SketchPatterns,SENNA-method]{SketchPatterns for SENNA}},
#'   \code{\link[=SketchPatterns,mSENNA-method]{SketchPatterns for mSENNA}}

setGeneric("SketchPatterns", function(senna, ...) {
  standardGeneric("SketchPatterns")
})



#' Sketch spatial patterns for selected SVGs
#'
#' Fit bin-wise weighted regression models for genes and extract their spatial trend patterns.  
#' This function processes each gene's expression profile along the curve axis, optionally within an interval.
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom tidyselect all_of
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom dplyr pick
#' @importFrom dplyr everything
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom methods new
#'
#' @param senna A `SENNA` object
#' @param genelist A character vector of genes to analyze. If `NULL`, all progression SVGs are included
#' @param bins Number of bins along the curve axis (optional if already binned)
#' @param active If `TRUE`, only active spots are used
#' @param interval Only values within this distance are included in regression
#' @param cores Number of cores for parallel processing
#'
#' @return An object of class `sktS` containing regression patterns and gene labels
#' @export

setMethod("SketchPatterns",
          signature = list(senna = "SENNA"),
          function(senna,
                   genelist = NULL,
                   bins = NULL,
                   active = FALSE,
                   interval = sqrt(2),
                   cores = 1L) {
            if(!"Bin" %in% colnames(senna@Coord[["Spatial"]])){
              senna <- GroupCurveCovariate(senna, bins = bins)
            }

            if(is.null(genelist)){
              genes <- c(senna@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]],
                        rev(senna@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]]))
            } else{
              rpt <- dplyr::filter(senna@Gene[["P.SVGs"]][["Report"]], Gene %in% genelist)
              pg <- split(rpt, rpt$Coefficients >= 0)[["TRUE"]]
              pg <- dplyr::arrange(pg, pg[["adj.p"]])[["Gene"]]
              ng <- split(rpt, rpt$Coefficients >= 0)[["FALSE"]]
              ng <- dplyr::arrange(ng, dplyr::desc(ng[["adj.p"]]))[["Gene"]]
              genes <- c(pg, ng)
            }



            gdf <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]),]
            if(active){
              idx <- rownames(gdf) %in% senna@Coord[["Activated"]][["cells"]]
              gdf <- gdf[idx, ]
            }
            gdf <- dplyr::select(gdf, tidyselect::all_of(genes))
            gdf <- dplyr::mutate(gdf, rn = rownames(gdf))
            cdf <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]),]
            if(active){
              cdf <- cdf[idx, ]
            }
            cdf <- dplyr::mutate(cdf, rn = rownames(cdf))

            oup <- parallel::mclapply(
              genes,
              function(svg) {
                inp <- merge(x = dplyr::select(gdf, tidyselect::all_of(svg) | rn),
                             y = dplyr::select(cdf, t | distance | Bin | rn),
                             by = "rn")
                inp <- dplyr::select(inp, -rn)
                inp <- dplyr::group_by(inp, Bin)
                res <- dplyr::summarise(inp, Coef = {
                  model <- stats::lm(get(svg) ~ t,
                                     weights = ifelse(distance <= interval, 1, 0),
                                     data = dplyr::pick(dplyr::everything()))
                  stats::coef(model)[[2]]
                })

                res <- dplyr::mutate(res, Gene = svg)
                return(res)
              }, mc.cores = cores)

            oup <- do.call(rbind, oup)

            oup <- methods::new(
              "sktS",
              sketch = oup,
              label = genes
            )

            return(oup)
          })



#' Sketch spatial patterns for selected SVGs (mSENNA)
#'
#' Fit bin-wise weighted regression models across multiple curve axes.  
#' Each gene's spatial profile is analyzed using a binned and weighted linear model.
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom tidyselect all_of
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom dplyr pick
#' @importFrom dplyr everything
#' @importFrom methods new
#'
#' @param senna A `mSENNA` object
#' @param genelist A character vector of genes to analyze. If `NULL`, all genes are included
#' @param msR A reference `msR` object, used to order genes by effect size (optional)
#' @param bins Number of bins to divide the curve axis (default = 4)
#' @param active If `TRUE`, only active spots are used
#' @param intervals A list of cutoff distances for each curve axis, used in weighted regression
#' @param cores Number of cores for parallel processing
#'
#' @return An object of class `sktS` containing regression patterns and gene labels
#' @export

setMethod("SketchPatterns",
          signature = list(senna = "mSENNA"),
          function(senna,
                   genelist = NULL,
                   msR = NULL,
                   bins = 4,
                   active = FALSE,
                   intervals,
                   cores = 1L) {

            slen <- length(senna@SENNA)
            vlen <- 1:slen
            if(is.null(genelist)){
              genes <- colnames(senna@SENNA[[1]]@Gene[["Spatial"]])
            } else{
              if(is.null(msR)) {
                genes <- genelist
                } else{
                  rpt <- dplyr::filter(msR@Report,
                                       Gene %in% genelist)
                  pg <- split(rpt, rpt$Coefficients >= 0)[["TRUE"]]
                  pg <- dplyr::arrange(pg, pg[["adj.p"]])[["Gene"]]
                  ng <- split(rpt, rpt$Coefficients >= 0)[["FALSE"]]
                  ng <- dplyr::arrange(ng, dplyr::desc(ng[["adj.p"]]))[["Gene"]]
                  genes <- c(pg, ng)
              }}


            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  t_scale = max(coord_mat$t) - min(coord_mat$t)
                  t_center = min(coord_mat$t)
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = (t - t_center) / t_scale,
                                             bat = v,
                                             interval = intervals[[v]])

                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  if(active){
                    idx <- rownames(gene_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    gene_mat <- gene_mat[idx, ]
                  }

                  return(gene_mat)
                }, mc.cores = cores)
            gmat <- do.call(rbind, gmat)

            gmat <- dplyr::mutate(gmat, rn = rownames(gmat))
            cmat <- dplyr::mutate(cmat, rn = rownames(cmat))
            Bin <- cut(cmat$t,
                       breaks = seq(0, 1, by = 1 / bins),
                       labels = FALSE,
                       include.lowest = TRUE)
            cmat$Bin <- factor(Bin)

            oup <- parallel::mclapply(
              genes,
              function(gene) {
                inp <- merge(x = dplyr::select(gmat, tidyselect::all_of(gene) | rn),
                             y = dplyr::select(cmat, t | distance | Bin | interval | rn),
                             by = "rn")
                inp <- dplyr::select(inp, -rn)
                inp <- dplyr::group_by(inp, Bin)
                res <- dplyr::summarise(inp, Coef = {
                  model <- stats::lm(get(gene) ~ t,
                                     weights = ifelse(distance <= interval, 1, 0),
                                     data = dplyr::pick(dplyr::everything()))
                  stats::coef(model)[[2]]
                })

                res <- dplyr::mutate(res, Gene = gene)
                return(res)
              }, mc.cores = cores)

            oup <- do.call(rbind, oup)

            oup <- methods::new(
              "sktS",
              sketch = oup,
              label = genes
            )

            return(oup)
          })



#' Plot a heatmap with SENNA-derived objects
#'
#' Define generic function `PatternPlot()` to visualize spatial gene expression patterns.  
#' This function dispatches to different methods depending on the input object type.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA`, `mSENNA`, `sktS`, or `simplesktS` object
#' @param ... Additional arguments passed to the method
#'
#' @seealso \code{\link[=PatternPlot,SENNA-method]{PatternPlot for SENNA}},   
#'   \code{\link[=PatternPlot,sktS-method]{PatternPlot for sktS}},  
#'   \code{\link[=PatternPlot,simplesktS-method]{PatternPlot for simplesktS}}
#'
#' @export

setGeneric("PatternPlot", function(senna, ...) {
  standardGeneric("PatternPlot")
})



#' Create a `simplesktS` object for SVG pattern sketching
#'
#' Generate a simple sketch matrix representing average expression patterns  
#' of selected genes along the curve axis, using either knot-centered or  
#' evenly binned segments depending on `nbins`.
#'
#' @importFrom methods new
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom stats complete.cases
#'
#' @param senna A `SENNA`, `msR`, or `sktS` object
#' @param genelist A character vector of gene names to include
#' @param nbins An integer specifying the number of bins for curve segmentation.  
#'   If `NULL`, the curve will be segmented around each knot using a fixed radius (`thr`).
#' @param interval The inclusion radius used to select spots near the curve axis (default = `sqrt(2)`)
#' @param thr Radius threshold used when `nbins = NULL` (default = `0.01`)
#' @param active If `TRUE`, restrict to activated spots as defined in `Coord$Activated`
#'
#' @return A `simplesktS` object containing the gene sketch matrix and spot count per segment
#' @export

Make_simplesktS <- function(senna,
                            genelist,
                            nbins = NULL,
                            interval = sqrt(2),
                            thr = 0.01,
                            active = FALSE){
  knt <- senna@CurveAxis[["knots"]]
  cord <- senna@Coord[["Spatial"]]
  gns <- dplyr::select(senna@Gene[["Spatial"]], 
                       all_of(genelist))
  
  if(is.null(nbins)){
    oup <- lapply(1:nrow(knt),
                  function(x){
                    nrm <- knt[x,]
                    cid <- (cord$X1 - nrm$X1)^2 +
                      (cord$X2 - nrm$X2)^2
                    cid <- (cid <= thr^2)
                    nsp <- sum(cid)
                    
                    mval <- colSums(gns[cid,]) / nsp
                    
                    return(c(mval, n_spots = nsp))
                  })
    
    oup <- do.call(rbind, oup)
    nsp <- as.integer(tibble::as_tibble(oup)$n_spots)
    oup <- dplyr::select(tibble::as_tibble(oup), genelist)
    oup <- t(oup)
    
    oup <- methods::new(
      "simplesktS",
      sketch = oup,
      label = genelist,
      n_spots = nsp
    )
  }
  
  else{
    ## Reflect curve length to CP.
    fx <- senna@CurveAxis[["fun"]][["x.coef"]]
    fy <- senna@CurveAxis[["fun"]][["y.coef"]]
    knots <- senna@CurveAxis[["knots"]]
    kn <- senna@CurveAxis[["fun"]][["t"]]
    
    if("spline" %in% senna@CurveAxis[["type"]]){
      if("trimmed" %in% senna@CurveAxis[["type"]]) {
        ll <- lapply(kn[-length(kn)],
                     function(k){
                       qf <- function(x){
                         sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                       }
                       return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                     })
        ll <- c(0, unlist(ll))
        cll <- cumsum(ll)
        tv <- cord[["t"]]
        tprime <- 
          ifelse(is.nan(tv),
                 tv,
                 {
                   fidx <- floor(tv)
                   fidx[fidx == 0] <- 1L
                   fidx[fidx > length(cll)] <- length(cll)
                   
                   (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                     cll[length(cll)]
                 })
        cord[["t"]] <- tprime
      } else {
        kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
        tv <- cord[["t"]]
        tmin <- min(tv); tmax <- max(tv)
        llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                        ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
        lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                        ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
        ll <- lapply(2:(mk-1),
                     function(k){
                       qf <- function(x){
                         sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                       }
                       return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                     })
        ll <- c(0, llw, unlist(ll), lhi)
        cll <- cumsum(ll)
        tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
        tv <- ifelse(tv > (max(kn) - 1L), 
                     (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                     tv)
        tv <- tv + 1L
        tprime <- {
          fidx <- floor(tv)
          fidx[fidx == 0] <- 1L
          fidx[fidx > length(cll)] <- length(cll)
          (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
            cll[length(cll)]
        }
        cord[["t"]] <- tprime
      }
    } else if("straight" %in% senna@CurveAxis[["type"]]){
      if("trimmed" %in% senna@CurveAxis[["type"]]) {
        ll <- lapply(kn[-length(kn)],
                     function(k){
                       dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                       dist <- sqrt(sum(dist))
                       return(dist)
                     })
        ll <- c(0, unlist(ll))
        cll <- cumsum(ll)
        tv <- cord[["t"]]
        tprime <- 
          ifelse(is.nan(tv),
                 tv,
                 {
                   fidx <- floor(tv)
                   fidx[fidx == 0] <- 1L
                   fidx[fidx > length(cll)] <- length(cll)
                   
                   (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                     cll[length(cll)]
                 })
        cord[["t"]] <- tprime
      } else {
        kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
        tv <- cord[["t"]]
        tmin <- min(tv); tmax <- max(tv)
        llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                        ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
        lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                        ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
        ll <- lapply(2:(mk-2),
                     function(k){
                       dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                       dist <- sqrt(sum(dist))
                       return(dist)
                     })
        ll <- c(0, llw, unlist(ll), lhi)
        cll <- cumsum(ll)
        tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
        tv <- ifelse(tv > (max(kn) - 1L), 
                     (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                     tv)
        tv <- tv + 1L
        tprime <- {
          fidx <- floor(tv)
          fidx[fidx == 0] <- 1L
          fidx[fidx > length(cll)] <- length(cll)
          
          (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
            cll[length(cll)]
        }
        cord[["t"]] <- tprime
      }
    }
    
    ## Evenly split
    gns <- gns[stats::complete.cases(cord),]
    cord <- cord[stats::complete.cases(cord),]
    cord$label <- cut(cord$t, 
                      breaks = seq(0, 1, length.out = nbins + 1), 
                      labels = 1:nbins, 
                      include.lowest = TRUE)
    gns <- gns[stats::complete.cases(cord),]
    cord <- cord[stats::complete.cases(cord),]
    
    oup <- lapply(1:nbins,
                  function(x){
                    idv <- with(cord, 
                                distance < interval &
                                  label == x)
                    if(active) idv <- idv * (rownames(cord) %in% senna@Coord[["Activated"]][["cells"]])
                    nsp <- sum(idv)
                    
                    mval <- colSums(gns[idv,]) / nsp
                    
                    return(c(mval, n_spots = nsp))
                  })
    
    oup <- do.call(rbind, oup)
    nsp <- as.integer(tibble::as_tibble(oup)$n_spots)
    oup <- dplyr::select(tibble::as_tibble(oup), genelist)
    oup <- t(oup)
    
    oup <- methods::new(
      "simplesktS",
      sketch = oup,
      label = genelist,
      n_spots = nsp
    )
    
  }}




#' Heatmap of gene patterns along the curve axis
#'
#' Plot a heatmap of gene expression patterns across curve segments  
#' using a `SENNA` object. Internally uses \code{\link{Make_simplesktS}}  
#' to generate a sketch matrix and delegates to  
#' \code{\link[=PatternPlot,simplesktS-method]{PatternPlot for simplesktS}}.
#'
#' @importFrom pheatmap pheatmap
#'
#' @param senna A `SENNA` object
#' @param genelist A character vector of genes to include in the heatmap.  
#'   If `NULL`, all genes in `Gene$Spatial` will be used.
#' @param n_bins Number of bins to divide the curve axis.  
#'   If `NULL`, expression will be averaged around each knot using `radius_thr`.
#' @param interval Curve-axis distance threshold to include spots (used only when `n_bins` is set)
#' @param radius_thr Radius for inclusion around knots (used only when `n_bins = NULL`)
#' @param xlab Logical. Show column (bin) labels. Default is `TRUE`.
#' @param ylab Logical. Show row (gene) labels. Default is `TRUE`.
#' @param legend Logical. Show legend in the heatmap. Default is `FALSE`.
#' @param colorset Optional color palette passed to `pheatmap::pheatmap()`
#'
#' @return A heatmap visualizing spatial gene expression patterns across curve segments
#' @export

setMethod("PatternPlot",
          signature = list(senna = "SENNA"),
          function(senna,
                   genelist = NULL,
                   n_bins = NULL,
                   interval = sqrt(2),
                   radius_thr = 0.01,
                   xlab = TRUE,
                   ylab = TRUE,
                   legend = FALSE,
                   colorset = NULL){
            if(is.null(genelist)){
              genes <- colnames(senna@Gene[["Spatial"]])
              }
            else {
              genes <- genelist
            }
            
            if(is.null(n_bins)){
              senna <- Make_simplesktS(senna,
                                       genelist = genes,
                                       thr = radius_thr)
            } else{
              senna <- Make_simplesktS(senna,
                                       genelist = genes,
                                       nbins = n_bins,
                                       interval = interval)
            }
            
            skt <- senna@sketch
            colnames(skt) <- seq(1, ncol(skt))

            if(is.null(colorset)){
              pheatmap::pheatmap(
                skt,
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                scale = "none",
                show_colnames = xlab,
                show_rownames = ylab,
                legend = legend)
            }
            
            else(
              pheatmap::pheatmap(
                skt,
                cluster_cols = FALSE,
                cluster_rows = TRUE,
                scale = "none",
                show_colnames = xlab,
                show_rownames = ylab,
                legend = legend,
                color = colorset)
            )
            
          })



#' Heatmap of spatial expression patterns
#'
#' Visualize gene expression patterns across curve segments using  
#' regression coefficients from a `sktS` object.  
#' The color scale highlights increasing and decreasing patterns.
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom scales rescale
#'
#' @param senna A `sktS` object containing spatial regression coefficients
#' @param colors A vector of five colors for gradient visualization,  
#'   typically from low to high (default: red–white–blue palette)
#'
#' @return A heatmap displaying increasing and decreasing expression trends for each gene
#' @export


setMethod("PatternPlot",
          signature = list(senna = "sktS"),
          function(senna,
                   colors = c("#8c1c1f", "#ffb6c1", "#fafafa", "#99d0e2", "#1c496e")){
            oup <- dplyr::mutate(senna@sketch,
                                 Gene = factor(Gene, levels = senna@label))
            ggplot2::ggplot(data = oup,
                            ggplot2::aes(Bin,
                                         Gene,
                                         fill = Coef)) +
              ggplot2::geom_raster() +
              ggplot2::scale_fill_gradientn(colors = colors,
                                            values = scales::rescale(c(min(oup$Coef),
                                                                       mean(oup$Coef[oup$Coef < 0]),
                                                                       0,
                                                                       mean(oup$Coef[oup$Coef > 0]),
                                                                       max(oup$Coef)),
                                                                     to = c(0, 1))) +
              ggplot2::theme(legend.position = "none",
                             axis.text.y = ggplot2::element_blank())
          }
)



#' Heatmap of mean expression near curve knots
#'
#' Visualize the average expression of selected genes  
#' near each curve knot using a `simplesktS` object.
#'
#' @importFrom pheatmap pheatmap
#'
#' @param senna A `simplesktS` object containing mean expression data
#' @param xlab Logical, whether to display x-axis labels (default: `TRUE`)
#' @param ylab Logical, whether to display y-axis labels (default: `TRUE`)
#' @param legend Logical, whether to display the color legend (default: `FALSE`)
#'
#' @return A heatmap showing the mean expression of selected genes near each knot
#' @export

setMethod("PatternPlot",
          signature = list(senna = "simplesktS"),
          function(senna,
                   xlab = TRUE,
                   ylab = TRUE,
                   legend = FALSE){
            skt <- senna@sketch
            colnames(skt) <- seq(1, ncol(skt))
            pheatmap::pheatmap(
              skt,
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              scale = "none",
              show_colnames = xlab,
              show_rownames = ylab,
              legend = legend
            )
          }
)


#' Show spatial bins along the curve axis
#'
#' Define generic function `ShowBins()` to visualize spatial bins  
#' segmented along the fitted curve axis.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` object
#' @param ... Additional arguments passed to the method
#'
#' @seealso \code{\link[=ShowBins,SENNA-method]{ShowBins for SENNA}}
#'
#' @export

setGeneric("ShowBins", function(senna, ...) {
  standardGeneric("ShowBins")
})





#' show spatial bins along the curve axis
#'
#' Visualize spatial bins assigned along a curve axis.  
#' The method divides the curve into evenly spaced bins (if specified)  
#' and shows how each spot is assigned based on its projected coordinate `t`.  
#' Spots outside the selected interval are drawn as background.
#'
#' @importFrom stats complete.cases
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_test
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna A `SENNA` object
#' @param bins Number of bins to divide the curve axis into
#' @param active Logical; whether to restrict to activated cells
#' @param interval Distance cutoff used to define proximity to the curve
#' @param colors Optional vector of colors to use for bin coloring
#' @param dot_size Size of dots for spots within interval
#' @param dot_alpha Transparency of dots within interval
#' @param bg_dot_color Color of background spots (outside interval)
#' @param bg_dot_size Size of background spots
#' @param bg_dot_alpha Transparency of background spots
#' @param line_color Color of curve axis
#' @param line_size Size of curve axis
#' @param line_alpha Transparency of curve axis
#' @param knots_color Color of knot points
#' @param knots_size Size of knot points
#' @param knots_alpha Transparency of knot points
#' @param order_label Logical; whether to annotate knot index numbers
#' @param text_size Size of the knot index labels
#'
#' @return A `ggplot` object showing spatial bin segmentation
#' @export

setMethod("ShowBins",
          signature = list(senna = "SENNA"),
          function(senna,
                   bins = NULL,
                   active = FALSE,
                   interval = sqrt(2),
                   colors = NULL,
                   dot_size = 1,
                   dot_alpha = 1,
                   bg_dot_color = "#cccccc",
                   bg_dot_size = 1,
                   bg_dot_alpha = 1,
                   line_color = "#7a7a7a",
                   line_size = 0.1,
                   line_alpha = 1,
                   knots_color = "#000000",
                   knots_size = 1.5,
                   knots_alpha = 1,
                   order_label = FALSE,
                   text_size = 3
          ){
            ## Reflect curve length to CP.
            cord <- senna@Coord[["Spatial"]]
            fx <- senna@CurveAxis[["fun"]][["x.coef"]]
            fy <- senna@CurveAxis[["fun"]][["y.coef"]]
            knots <- senna@CurveAxis[["knots"]]
            kn <- senna@CurveAxis[["fun"]][["t"]]
            
            if("spline" %in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- cord[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
                cord[["t"]] <- tprime
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- cord[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                ll <- lapply(2:(mk-1),
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx > length(cll)] <- length(cll)
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
                cord[["t"]] <- tprime
              }
            } else if("straight" %in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- cord[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
                cord[["t"]] <- tprime
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- cord[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                ll <- lapply(2:(mk-2),
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx > length(cll)] <- length(cll)
                  
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
                cord[["t"]] <- tprime
              }
            }
            
            ## Evenly split
            cord <- cord[stats::complete.cases(cord),]
            cord$bin <- cut(cord$t, 
                            breaks = seq(0, 1, length.out = bins + 1), 
                            labels = 1:bins, 
                            include.lowest = TRUE)
            cord <- cord[stats::complete.cases(cord),]
            cord$bin <- factor(cord$bin)
            cord <- cord[cord$distance <= interval,]
            
            bg <- setdiff(
              senna@Coord[["Spatial"]][, c("X1", "X2")], 
              cord[, c("X1", "X2")])
            
            p <- ggplot2::ggplot() +
              ggplot2::geom_point(
                ggplot2::aes(X1, X2), data = bg,
                color = bg_dot_color, size = bg_dot_size, alpha = bg_dot_alpha) +
              ggplot2::geom_point(
                ggplot2::aes(X1, X2, color = bin), data = cord,
                size = dot_size, alpha = dot_alpha)
            
            if(!is.null(colors)) {
              p <- p +
                ggplot2::scale_color_manual(values = colors)
            }
            
            p <- p +
              ggplot2::geom_point(
                ggplot2::aes(X1, X2), data = knots,
                color = knots_color, size = knots_size, alpha = knots_alpha) +
              ggplot2::geom_point(
                ggplot2::aes(X1, X2), data = crvtrjry(senna)) +
              ggplot2::theme_test()
            
            if(order_label){
              p <- p +
                ggrepel::geom_text_repel(
                  ggplot2::aes(X1, X2),
                  label = 1:nrow(knots),
                  data = knots,
                  max.overlaps = nrow(knots),
                  size = text_size)
            }
            
            p
          })





#' Show spatial distribution within a radius threshold
#'
#' Visualize which spots fall within a radius threshold around each knot  
#' on the curve axis. Spots within the threshold are grouped by knot index  
#' and optionally colored. Remaining spots are shown as background.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_test
#' @importFrom ggplot2 labs
#' @importFrom dplyr mutate
#'
#' @param senna A `SENNA` object
#' @param radius_thr A numeric radius threshold used to define proximity
#' @param inner_color Optional vector of colors for groups within the threshold
#' @param outer_color Color for spots outside the radius (default = "#bfbfbf")
#' @param dots_alpha Transparency of all spots
#' @param dots_size Size of all spots
#' @param line_color Color of the curve axis
#' @param line_size Size of the curve axis
#' @param line_alpha Transparency of the curve axis
#' @param knots_color Color of the knot points
#' @param knots_size Size of the knot points
#' @param knots_alpha Transparency of the knot points
#'
#' @return A `ggplot` object visualizing radius-threshold-based segmentation
#' @export

ShowRadiusThr <- function(senna,
                          radius_thr,
                          inner_color = NULL,
                          outer_color = "#bfbfbf",
                          dots_alpha = 1,
                          dots_size = 1,
                          line_color = "#7a7a7a",
                          line_size = 0.2,
                          line_alpha = 1,
                          knots_color = "#000000",
                          knots_size = 2,
                          knots_alpha = 1){
  knt <- senna@CurveAxis[["knots"]]
  cord <- senna@Coord[["Spatial"]][, 1:2]
  griddat <- crvtrjry(senna)

  ids <- lapply(1:nrow(knt),
                function(k){
                  nrm <- knt[k,]
                  oup <- dplyr::mutate(cord,
                                       D = (cord$X1 - nrm$X1)**2 + (cord$X2 - nrm$X2)**2,
                                       K = k)
                  oup <- oup[oup$D <= radius_thr**2,]
                  return(oup)
                })

  ids <- do.call(rbind, ids)
  ods <- cord[!rownames(cord) %in% rownames(ids),]

  if(is.null(inner_color)){
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(X1, X2, col = as.factor(K)), data = ids,
                          size = dots_size, alpha = dots_alpha)
  } else{
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(X1, X2, col = as.factor(K)), data = ids,
                          size = dots_size, alpha = dots_alpha) +
      ggplot2::scale_color_manual(values = inner_color)
  }


  p <- p +
    ggplot2::geom_point(ggplot2::aes(X1, X2), data = ods,
                        col = outer_color, size = dots_size, alpha = dots_alpha) +
    ggplot2::geom_point(ggplot2::aes(X1, X2), data = griddat,
                        col = line_color, size = line_size, alpha = line_alpha) +
    ggplot2::geom_point(ggplot2::aes(X1, X2), data = knt,
                        col = knots_color, size = knots_size, alpha = knots_alpha) +
    ggplot2::theme_test() +
    ggplot2::labs(color = "Order")

  p
}



#' Volcano plot for peak-based SVGs
#'
#' Define generic function `PeakVolPlot()`.
#'
#' @importFrom methods setGeneric
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments passed to methods
#' @return A volcano plot or list of plots highlighting peak-based SVGs
#' @export
#' @seealso \code{\link{PeakVolPlot,SENNA-method}}

setGeneric("PeakVolPlot", function(senna, ...) {
  standardGeneric("PeakVolPlot")
})



#' show volcano plots for peak-based SVGs
#'
#' Visualize progression or regionation SVGs identified at specific peaks  
#' along the curve axis. Separate volcano plots are created for each peak.
#'
#' @importFrom methods setMethod
#' @importFrom dplyr arrange
#' @importFrom dplyr setdiff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggrepel geom_text_repel
#'
#' @param senna A `SENNA` object containing peak-based SVGs
#' @param FDR_level FDR threshold used to define significance
#' @param grad_cutoff Threshold for gradient coefficient magnitude
#' @param dot_size Size of plot dots
#' @param dot_alpha Transparency of plot dots
#' @param positive Color for positively associated genes
#' @param negative Color for negatively associated genes
#' @param cutoff_col Color for significance threshold lines
#' @param nrepel Number of top genes to label with text
#' @param text_size Size of text labels
#' @param yfix Whether to fix y-axis limit
#'
#' @return A list of `ggplot` objects, one for each peak
#' @exportMethod PeakVolPlot

setMethod("PeakVolPlot",
          signature(senna = "SENNA"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0,
                   dot_size = 0.5,
                   dot_alpha = 0.5,
                   positive = "#ec1b24",
                   negative = "#1b74bc",
                   cutoff_col = "#888888",
                   nrepel = 0L,
                   text_size = 3,
                   yfix = FALSE){

            if(senna@Gene[["PeakSVGs"]][["type"]] == "Progression"){
              reports <- senna@Gene[["PeakSVGs"]][["Reports"]]
              len <- length(reports)
              ps <- lapply(
                1:len,
                function(n){
                  k <- reports[[n]][["Peak"]]
                  report <- reports[[n]][["Report"]]
                  possig <- subset(report, Gene %in% reports[[n]][["Variable_gene"]][["positive"]])
                  negsig <- subset(report, Gene %in% reports[[n]][["Variable_gene"]][["negative"]])
                  nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

                  if(yfix){
                    p <- ggplot2::ggplot() +
                      ggplot2::ylim(c(NA, 301)) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = 'gray', alpha = dot_alpha,
                                          size = dot_size, data = nonsig) +
                      ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                           col = positive, alpha = dot_alpha,
                                           size = dot_size, data = possig) +
                      ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                           col = negative, alpha = dot_alpha,
                                           size = dot_size, data = negsig) +
                      ggplot2::labs(x = "Gradients (Scaled)",
                                    y = "Adjusted p-value, -log",
                                    title = paste0("Peak = ", k)) +
                      ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                      ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                      ggplot2::theme_test() +
                      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())
                  } else{
                    p <- ggplot2::ggplot() +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = 'gray', alpha = dot_alpha,
                                          size = dot_size, data = nonsig) +
                      ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                           col = positive, alpha = dot_alpha,
                                           size = dot_size, data = possig) +
                      ggplot2::geom_jitter(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                           col = negative, alpha = dot_alpha,
                                           size = dot_size, data = negsig) +
                      ggplot2::labs(x = "Gradients (Scaled)",
                                    y = "Adjusted p-value, -log",
                                    title = paste0("Peak = ", k)) +
                      ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                      ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                      ggplot2::theme_test() +
                      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())
                  }

                  if(grad_cutoff != 0){
                    p <- p +
                      ggplot2::geom_vline(xintercept = sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col) +
                      ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff), linetype = "dotted", alpha = 0.7, linewidth = 1.3, col = cutoff_col)
                  }

                  if(nrepel > 0){
                    ntop <- dplyr::arrange(report, log(adj.p))[1:nrepel, ]
                    p <- p +
                      ggrepel::geom_text_repel(
                        aes(Coefficients / sqrt(abs(Coefficients)),
                            -log10(adj.p + 1e-300),
                            label = Gene),
                        size = text_size,
                        data = ntop)
                  }

                  return(p)
                  })
            }

            else if(senna@Gene[["PeakSVGs"]][["type"]] == "Regionation"){
              reports <- senna@Gene[["PeakSVGs"]][["Reports"]]
              len <- length(reports)
              ps <- lapply(
                1:len,
                function(n){
                  k <- reports[[n]][["Peak"]]
                  report <- reports[[n]][["Report"]]
                  possig <- subset(report, Gene %in% reports[[n]][["Variable_gene"]][["positive"]])
                  negsig <- subset(report, Gene %in% reports[[n]][["Variable_gene"]][["negative"]])
                  nonsig <- dplyr::setdiff(report, rbind(possig, negsig))

                  if(yfix){
                    p <- ggplot2::ggplot() +
                      ggplot2::ylim(c(NA, 301)) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = 'gray', alpha = dot_alpha,
                                          size = dot_size, data = nonsig) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = positive, alpha = dot_alpha,
                                          size = dot_size, data = possig) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = negative, alpha = dot_alpha,
                                          size = dot_size, data = negsig) +
                      ggplot2::labs(x = "Gradients (Scaled)",
                                    y = "Adjusted p-value, -log",
                                    title = paste0("Peak = ", k)) +
                      ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                      ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                      ggplot2::theme_light() +
                      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())
                  } else{
                    p <- ggplot2::ggplot() +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = 'gray', alpha = dot_alpha,
                                          size = dot_size, data = nonsig) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = positive, alpha = dot_alpha,
                                          size = dot_size, data = possig) +
                      ggplot2::geom_point(aes(Coefficients / sqrt(abs(Coefficients)), -log10(adj.p + 1e-300)),
                                          col = negative, alpha = dot_alpha,
                                          size = dot_size, data = negsig) +
                      ggplot2::labs(x = "Gradients (Scaled)",
                                    y = "Adjusted p-value, -log",
                                    title = paste0("Peak = ", k)) +
                      ggplot2::geom_vline(xintercept = 0, linetype = "dotdash", alpha = 0.7) +
                      ggplot2::geom_hline(yintercept = -log10(FDR_level + 1e-300), linetype = "dashed", alpha = 0.7, col = cutoff_col, linewidth = 1.3) +
                      ggplot2::theme_light() +
                      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank())
                  }

                  if(grad_cutoff != 0){
                    p <- p +
                      ggplot2::geom_vline(xintercept = sqrt(grad_cutoff),
                                          linetype = "dotted", alpha = 0.7,
                                          linewidth = 1.3, col = cutoff_col) +
                      ggplot2::geom_vline(xintercept = -sqrt(grad_cutoff),
                                          linetype = "dotted", alpha = 0.7,
                                          linewidth = 1.3, col = cutoff_col)
                  }

                  if(nrepel > 0){
                    ntop <- dplyr::arrange(report, log(adj.p))[1:nrepel, ]
                    p <- p +
                      ggrepel::geom_text_repel(
                        aes(Coefficients / sqrt(abs(Coefficients)),
                            -log10(adj.p + 1e-300),
                            label = Gene),
                        size = text_size,
                        data = ntop)
                  }})

              return(ps)
            }})



