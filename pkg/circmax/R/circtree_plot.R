plot.circtree <- function(x, terminal_panel = node_circular,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...){

  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    partykit:::plot.constparty(partykit::as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    partykit::plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}


plot_circular <- function(X, coefs = NULL, stack = 2, cex = 0.8, label = TRUE, 
  circlab = c('$0$', '$\\pi/2$', '$\\pi$','$3/2\\pi$'), fit_dens = FALSE){

  par(mar = c(0, 0, 0, 0) + 0., cex = cex)

  circ_ln <- seq(0, 360, 0.5) * pi / 180
  plot(cos(circ_ln), sin(circ_ln), xlim = c(-2,2), ylim = c(-2,2),
    type = "l", axes = FALSE, xlab = "",ylab = "", asp = 1)
  
  circ_sgm <- seq(0, 360, 45) * pi / 180
  segments(.9 * cos(circ_sgm), .9 * sin(circ_sgm), 0.99 * cos(circ_sgm), 0.99 * sin(circ_sgm))

  if(label & requireNamespace("latex2exp", quietly = TRUE)){
    text(.7, 0, latex2exp::TeX(circlab[1]))
    text(0, .7, latex2exp::TeX(circlab[2]))
    text(-.7, 0, latex2exp::TeX(circlab[3]))
    text(0, -.7, latex2exp::TeX(circlab[4]))
  }

  Xhist <- hist(X, plot = FALSE, breaks = seq(0, 360, stack) * pi / 180)

  idx <- (Xhist$density != 0)
  points(1 * cos(Xhist$mids), 1 * sin(Xhist$mids), pch = 20, col = c(NA, gray(0.2, alpha = 0.7))[idx + 1]) 
  
  Xscaled <- Xhist$density / max(Xhist$density)
  for(i in seq_along(Xhist$mids)) {
    lines(c(1 * cos(Xhist$mids[i]), (Xscaled[i] + 1) * cos(Xhist$mids[i])), 
      c(1 * sin(Xhist$mids[i]), (Xscaled[i] + 1) * sin(Xhist$mids[i])), col = gray(0.5))
  }
  
  if(fit_dens){  
    d <- circular::density.circular(circular::circular(X), bw = 1)
    lines(d, col = 2, lty = 2, lwd = 1.5)
  }
  if(!is.null(coefs)){
    arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0.1, code = 2, lwd = 2)
    circular::plot.function.circular(function(x) 
      dvonmises(x, circular::circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5)
  }
}


node_circular <- function(obj, which = NULL, id = TRUE, pop = TRUE,
  xlab = FALSE, ylab = FALSE, mainlab = NULL, ...){

  ## obtain dependent variable
  y <- obj$fitted[["(response)"]]
  fitted <- obj$fitted[["(fitted)"]]

  ## y- and x-lab
  if(isTRUE(ylab)) ylab <- names(y)
  if(identical(ylab, FALSE)) ylab <- ""
  if(identical(xlab, FALSE)) xlab <- ""

  ## set up appropriate panel functions
  rval <- function(node) {
    
    ## node index
    nid <- partykit::id_node(node)
    ix <- fitted %in% partykit::nodeids(obj, from = nid, terminal = TRUE)

    ## dependent variable
    y <- y[ix]
    if(length(y) > 1000) y <- sample(y, 1000)

    ## coefficients
    coefs <- c(signif(coef(node$info$object)[1],2), signif(coef(node$info$object)[2],2))

    ## set up viewport tree
    top.vp <- grid::viewport(layout = grid::grid.layout(ncol = 1, nrow = 3,
		       widths = grid::unit(c(1, 1, 1), c("null", "null", "null")), heights = grid::unit(c(0.3, 0.1, 0.6), c("npc", "npc", "npc"))),
		       width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), name = paste("node", nid, sep = ""))
    topmargin <- grid::viewport(layout.pos.col = 1, layout.pos.row = 1, name = "topmargin")
    text <- grid::viewport(layout.pos.col = 1, layout.pos.row = 2, name = "text")
    plot <- grid::viewport(layout.pos.col = 1, layout.pos.row = 3, name = "plot", )
    
    splot <- grid::vpTree(top.vp, grid::vpList(topmargin, text, plot))
    grid::pushViewport(splot)

    ## main title
    grid::seekViewport("text")

    if (is.null(mainlab)) { 
      mainlab <- if(id) {
    	function(id, nobs, mu, kappa) sprintf("Node %s (n = %s) \n mu = %s, kappa = %s", id, nobs, mu, kappa)
    	#function(id, nobs, mu, kappa) sprintf("Node %s (n = %s)", id, nobs)
      } else {
    	function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, node$info$object$ny, coefs[1], coefs[2])
    }
    grid::grid.text(mainlab, y = grid::unit(0.9, "npc"))

    ## plot rectangle and actual graphic, optional x and ylab 
    grid::seekViewport("plot")

    grid::grid.rect(gp = grid::gpar(fill = "transparent", col = 1), width = grid::unit(0.9, "npc"))
    gridGraphics::grid.echo(function() plot_circular(y, coefs, ...), newpage = FALSE)

    if(ylab != "") grid::grid.text(ylab, y = grid::unit(0.5, "npc"), x = grid::unit(-2.5, "lines"), rot = 90)
    if(xlab != "") grid::grid.text(xlab, x = grid::unit(0.5, "npc"), y = grid::unit(-2, "lines")) 	       

    if(pop) grid::popViewport() else grid::upViewport()
    if(pop) grid::popViewport() else grid::upViewport()
  }
  
  return(rval)
}
class(node_circular) <- "grapcon_generator"

