circtree_plot <- function(x, terminal_panel = node_circular,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot.constparty(as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    partykit::plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}


# Small helper functions
plot_circ <- function(x, plt_stack = TRUE, plt_rose = FALSE, coefs = NULL, rotation = "clock", zero = pi/2, ...){
  require("circular", quietly = TRUE)
  suppressWarnings(x <- as.circular(x, rotation = rotation, zero = zero))
  if(plt_stack){
    plot(x, cex = 0.8, control.circle = circle.control(cex = 0.5), bin = 720, stack = TRUE, sep = 0.035, 
      rotation = rotation, zero = zero, col = gray(0.25, alpha = 0.2), ...)
    #axis.circular(at = circular(seq(0, 7*pi/4, pi/2)), labels = c("N", "W", "E", "S"))
  } else{ 
    rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, rotation = rotation, zero = zero, axes = FALSE, ...)
  }
  if(plt_rose){
    rose.diag(x, bins = 16, col = "darkgrey", cex = 1., prop = 1.3, add = TRUE, rotation = rotation, zero = zero, axes = FALSE)
  }
  #lines(density.circular(x, bw = coefs[2]), lty = 2)
  if(!is.null(coefs)){
    circular::arrows.circular(circular(coefs[1]), 0.4, col = 2, lwd = 2, length = 0.05, rotation = rotation, zero = zero, code = 2)
    plot.function.circular(function(x) dvonmises(x, circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5,
      rotation = rotation, zero = zero)
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
    top_vp <- viewport(layout = grid.layout(nrow = 1, ncol = 2,
		       widths = unit(c(ylines, 1), c("lines", "null")), heights = unit(1, "null")),
		       width = unit(1, "npc"), height = unit(1, "npc")  - unit(2, "lines"),
		       name = paste("node", nid, sep = ""))
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))

    ## main title
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)

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
    grid.text(mainlab, y = unit(1, "npc") - unit(1.2, "lines"))
    popViewport()

    ## select panel
    plot_vpi <- viewport(layout.pos.col = 2L, layout.pos.row = 1L)
    pushViewport(plot_vpi)

    ## call panel function
    name <- paste("node", nid, "-", sep = "")
    pushViewport(plotViewport(margins = margins, name = name))
    grid.rect(gp = gpar(fill = "transparent", col = 1))
    gridGraphics::grid.echo(function() plot_circ(y, plt_stack = plt_stack, plt_rose = plt_rose, coefs = coefs, 
      rotation = rotation, zero = zero, shrink = shrink, tol = tol, tcl.text = tcl.text, ...), newpage = FALSE)
    #gridGraphics::grid.echo(function() plot(y), newpage = FALSE)

    #sublab <- function(mu, kappa) sprintf("mu = %s, kappa = %s", mu, kappa)
    #grid.text(sublab(coefs[1], coefs[2]), y = unit(0, "npc"), x = unit(0.5, "npc"))
    if(ylab != "") grid.text(ylab, y = unit(0.5, "npc"), x = unit(-2.5, "lines"), rot = 90)
    if(xlab != "") grid.text(xlab, x = unit(0.5, "npc"), y = unit(-2, "lines")) 	       

    if(pop) popViewport() else upViewport()
    if(pop) popViewport() else upViewport()
    if(pop) popViewport() else upViewport()
  }
  
  return(rval)
}
class(node_circular) <- "grapcon_generator"