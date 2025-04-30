#' A function to draw qq plot of a z-residual
#'
#' @param legend.args args
#' @export

legend.box <- function(legend.args){
  # Get legend dimensions without plotting it
  legend_info <- do.call(legend, c(legend.args, list(plot = FALSE)))

  args.new <- modifyList(legend.args, list(bty = "n"))
  do.call(legend, args.new)

  legend.text <- args.new[["legend"]]
  legend.cex <- args.new[["cex"]]
  # Draw a custom box around the legend using the extracted dimensions
  # rect(xleft = legend_info$rect$left,
  #      ybottom = legend_info$rect$top - legend_info$rect$h,
  #      xright = legend_info$rect$left + legend_info$rect$w,
  #      ytop = legend_info$rect$top,
  #      border = "black", lwd = 0.5)

  rect(xleft = legend_info$text$x[1]-1,
       ybottom = legend_info$text$y[length(legend_info$text$y)]-legend.cex,
       xright = legend_info$text$x[length(legend_info$text$x)] + strwidth(legend.text, cex = legend.cex)[length(legend.text)],
       ytop = legend_info$text$y[1]+legend.cex,
       border = "black", lwd = 0.5)
}
