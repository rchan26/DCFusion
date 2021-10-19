#' @export
obtain_hypercube_vertices <- function(bessel_layers) {
  if (!is.list(bessel_layers)) {
    stop("hypercube_vertices: bessel_layers must be a list")
  }
  bounds <- lapply(1:length(bessel_layers), function(d) c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
  B <- lapply(1:length(bessel_layers), function(d) {
    if ((bessel_layers[[d]]$L < 0) & (bessel_layers[[d]]$U > 0)) {
      return(c(bessel_layers[[d]]$L, 0, bessel_layers[[d]]$U))
    } else {
      return(c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
    }
  })
  vertices <- as.matrix(expand.grid(bounds))
  V <- as.matrix(expand.grid(B))
  colnames(vertices) <- c()
  colnames(V)
  return(list('vertices' = vertices, 'V' = V))
}

#' #' @export
#' surrounding_hypercube_vertices <- function(hypercube_vertices) {
#'   bounds <- lapply(1:ncol(hypercube_vertices), function(d) c(min(hypercube_vertices[,d]), max(hypercube_vertices[,d])))
#'   vertices <- as.matrix(expand.grid(bounds))
#'   colnames(vertices) <- c()
#'   return(vertices)
#' }
