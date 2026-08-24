#' A dashed reference line that's hoverable along its entire length, not just its endpoints
#'
#' `ggplot2::geom_vline()`/`geom_hline()` are drawn as a single 2-point segment, so
#' `plotly::ggplotly()`'s hover (which snaps to the nearest data point on a line trace) only
#' triggers near the two endpoints -- everywhere else along the line, nothing happens. Using many
#' points along the same segment instead means there's always a nearby point to snap to, so hover
#' works no matter where along the line the cursor is.
#'
#' @param p a `ggplot2` object, with its data/scale layers already added (this function reads
#'   their actual rendered range to span the line fully, so call it last).
#' @param orientation `"vertical"` (constant x, spanning the panel's y range) or `"horizontal"`
#'   (constant y, spanning the panel's x range).
#' @param at the x (if vertical) or y (if horizontal) position of the line.
#' @param text hover text for the line.
#' @param n number of points along the line (100 is dense enough that no gap is perceptible).
#' @return `p` with the reference line layer added.
#' @keywords internal
hover_line <- function(p, orientation = c("vertical", "horizontal"), at, text, n = 100) {
  orientation <- match.arg(orientation)
  built <- ggplot2::ggplot_build(p)
  xrange <- built$layout$panel_scales_x[[1]]$range$range
  yrange <- built$layout$panel_scales_y[[1]]$range$range
  line_df <- if (orientation == "vertical") {
    data.frame(x = at, y = seq(yrange[1], yrange[2], length.out = n), text = text)
  } else {
    data.frame(x = seq(xrange[1], xrange[2], length.out = n), y = at, text = text)
  }
  # ggplot2 warns "Ignoring unknown aesthetics: text" since geom_line() has no
  # such aesthetic natively -- harmless; the value still reaches the built
  # layer data that plotly::ggplotly() reads it from (verified empirically).
  p + suppressWarnings(ggplot2::geom_line(
    data = line_df, mapping = ggplot2::aes(x = x, y = y, text = text),
    inherit.aes = FALSE, linetype = "dashed", color = "red"
  ))
}
