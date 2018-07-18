plot.circtree <- function(x, terminal_panel = node_circular,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
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


## Small helper functions
#plot_circ <- function(x, plt_stack = TRUE, plt_rose = FALSE, coefs = NULL, rotation = "clock", zero = pi/2, ...){
#  suppressWarnings(x <- circular::as.circular(x, rotation = rotation, zero = zero))
#  if(plt_stack){
#    plot(x, cex = 0.8, control.circle = circular::circle.control(cex = 0.5), bin = 720, stack = TRUE, sep = 0.035, 
#      rotation = rotation, zero = zero, col = gray(0.25, alpha = 0.2), ...)
#    #axis.circular(at = circular(seq(0, 7*pi/4, pi/2)), labels = c("N", "W", "E", "S"))
#  } else{ 
#    circular::rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, rotation = rotation, zero = zero, axes = FALSE, ...)
#  }
#  if(plt_rose){
#    circular::rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, add = TRUE, rotation = rotation, zero = zero, axes = FALSE)
#  }
#  #lines(density.circular(x, bw = coefs[2]), lty = 2)
#  if(!is.null(coefs)){
#    circular::arrows.circular(circular::circular(coefs[1]), 0.4, col = 2, lwd = 2, length = 0.05, rotation = rotation, zero = zero, code = 2)
#    circular::plot.function.circular(function(x) dvonmises(x, circular::circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5,
#      rotation = rotation, zero = zero)
#  }
#}


plot_circular <- function(X, coefs = NULL, stack = 2, cex = 0.8, label = TRUE){

  par(mar = c(0, 0, 2, 0) + 0.1, cex = cex)

  circ_ln <- seq(0, 360, 0.5) * pi / 180
  plot(cos(circ_ln), sin(circ_ln), xlim = c(-2,2), ylim = c(-2,2),
    type = "l", axes = FALSE, xlab = "",ylab = "", asp = 1)
  
  circ_sgm <- seq(0, 360, 45) * pi / 180
  segments(.9 * cos(circ_sgm), .9 * sin(circ_sgm), 0.99 * cos(circ_sgm), 0.99 * sin(circ_sgm))

  if(label & requireNamespace("latex2exp", quietly = TRUE)){
    #text(.7, 0, latex2exp::TeX('$0 / 2\\pi$'))
    #text(-.7, 0, expression(pi))
    #text(0, -.7, latex2exp::TeX('$\\frac{3\\,\\pi}{2}$'))
    #text(0, .7, latex2exp::TeX('$\\frac{\\pi}{2}$'))
    text(.7, 0, latex2exp::TeX('$0$'))
    text(-.7, 0, expression(pi))
    text(0, -.7, latex2exp::TeX('$3/2\\,\\pi$'))
    text(0, .7, latex2exp::TeX('$\\pi/2$'))
  }

  Xhist <- hist(X, plot = FALSE, breaks = seq(0, 360, stack) * pi / 180)

  idx <- (Xhist$density != 0)
  points(1 * cos(Xhist$mids), 1 * sin(Xhist$mids), pch = 20, col = c(NA, gray(0.2, alpha = 0.7))[idx + 1]) 
  
  Xscaled <- Xhist$density / max(Xhist$density)
  for(i in seq_along(Xhist$mids)) {
    lines(c(1 * cos(Xhist$mids[i]), (Xscaled[i] + 1) * cos(Xhist$mids[i])), 
      c(1 * sin(Xhist$mids[i]), (Xscaled[i] + 1) * sin(Xhist$mids[i])), col = gray(0.5))
  }
  
  #d <- circular::density.circular(X, bw = 2)
  #lines(d, col = 2)
  if(!is.null(coefs)){
    arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0.1, code = 2, lwd = 2)
    circular::plot.function.circular(function(x) 
      dvonmises(x, circular::circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5)
  }
  
}


node_circular <- function(obj, which = NULL, id = TRUE, pop = TRUE,
  bg = "white", ylines = NULL, xlab = FALSE, ylab = FALSE, margins = rep(0.2, 4),
  mainlab = NULL, rotation = "clock", zero = pi/2, plt_stack = TRUE, plt_rose = FALSE, shrink = 1, tol = 0.04, tcl.text = 0.125, ...)
  {

  ## obtain dependent variable
  y <- obj$fitted[["(response)"]]
  fitted <- obj$fitted[["(fitted)"]]

  ## y- and x-lab
  if(isTRUE(ylab)) ylab <- names(y)
  if(identical(ylab, FALSE)) ylab <- ""
  if(is.null(ylines)) ylines <- ifelse(identical(ylab, ""), 0, 2)
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

    ## set up top viewport
    top_vp <- grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2,
		       widths = grid::unit(c(ylines, 1), c("lines", "null")), heights = grid::unit(1, "null")),
		       width = grid::unit(1, "npc"), height = grid::unit(1, "npc")  - grid::unit(4, "lines"),
		       name = paste("node", nid, sep = ""))
    grid::pushViewport(top_vp)
    grid::grid.rect(gp = grid::gpar(fill = bg, col = 0))

    ## main title
    top <- grid::viewport(layout.pos.col = 2, layout.pos.row = 1)
    grid::pushViewport(top)

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
    grid::grid.text(mainlab, y = grid::unit(1, "npc") - grid::unit(1.4, "lines"))
    grid::popViewport()

    ## select panel
    plot_vpi <- grid::viewport(layout.pos.col = 2L, layout.pos.row = 1L)
    grid::pushViewport(plot_vpi)

    ## call panel function
    name <- paste("node", nid, "-", sep = "")
    grid::pushViewport(grid::plotViewport(margins = margins, name = name, xscale = c(-2, 2), yscale = c(-2, 2)))
    grid::grid.rect(gp = grid::gpar(fill = "transparent", col = 1))
    #gridGraphics::grid.echo(function() plot_circ(y, plt_stack = plt_stack, plt_rose = plt_rose, coefs = coefs, 
    #  rotation = rotation, zero = zero, shrink = shrink, tol = tol, tcl.text = tcl.text, ...), newpage = FALSE)
    #gridGraphics::grid.echo(function() plot(y), newpage = FALSE)
    gridGraphics::grid.echo(function() plot_circular(y, coefs), newpage = FALSE)

    #sublab <- function(mu, kappa) sprintf("mu = %s, kappa = %s", mu, kappa)
    #grid::grid.text(sublab(coefs[1], coefs[2]), y = grid::unit(0, "npc"), x = grid::unit(0.5, "npc"))
    if(ylab != "") grid::grid.text(ylab, y = grid::unit(0.5, "npc"), x = grid::unit(-2.5, "lines"), rot = 90)
    if(xlab != "") grid::grid.text(xlab, x = grid::unit(0.5, "npc"), y = grid::unit(-2, "lines")) 	       

    if(pop) grid::popViewport() else grid::upViewport()
    if(pop) grid::popViewport() else grid::upViewport()
    if(pop) grid::popViewport() else grid::upViewport()
  }
  
  return(rval)
}
class(node_circular) <- "grapcon_generator"

