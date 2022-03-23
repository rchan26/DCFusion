#' Hypercube vertices
#'
#' Obtains the vertices of the hypercube given Bessel layer information
#'
#' @param bessel_layers a list of length dim where list[[d]] is the Bessel layer
#'                      for component d, which is a list of length 4
#'                      where L denotes the hard lower bound, l denotes the soft
#'                      lower bound, u denotes the soft upper bound and U 
#'                      denotes the hard upper bound
#' @param dim dimension
#'
#' @return A list of length 2 with
#' \describe{
#'  \item{vertices}{matrix where each row is a vertex of the hypercube}
#'  \item{V}{matrix where each row is a vertex of the hypercube including points
#'           crossing the origin in each dimension}
#' }
#'
#' @export
obtain_hypercube_vertices <- function(bessel_layers, dim) {
  if (!is.list(bessel_layers)) {
    stop("hypercube_vertices: bessel_layers must be a list of length dim")
  }
  if (dim == 1) {
    if (!identical(names(bessel_layers), c("L", "l", "u", "U"))) {
      stop("hypercube_vertices: if dim==1, bessel_layers must be a list of length 4 with names (L, l, u, U)")
    }
    return(list('vertices' = matrix(c(bessel_layers$L, bessel_layers$U)),
                'V' = matrix(c(bessel_layers$L, bessel_layers$U))))
  } else if (dim > 1) {
    if (length(bessel_layers)!=dim) {
      stop("hypercube_vertices: if dim > 1, bessel_layers must be a list of length dim")
    } else if (!all(sapply(1:dim, function(d) identical(names(bessel_layers[[d]]), c("L", "l", "u", "U"))))) {
      stop("hypercube_vertices: if dim > 1, bessel_layers[[d]] must be a list of length 4 with names (L, l, u, U)")
    }
    bounds <- lapply(1:dim, function(d) c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
    B <- lapply(1:dim, function(d) {
      if ((bessel_layers[[d]]$L < 0) & (bessel_layers[[d]]$U > 0)) {
        return(c(bessel_layers[[d]]$L, 0, bessel_layers[[d]]$U))
      } else {
        return(c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
      }
    })
    vertices <- as.matrix(expand.grid(bounds))
    V <- as.matrix(expand.grid(B))
    colnames(vertices) <- c()
    colnames(V) <- c()
    return(list('vertices' = vertices, 'V' = V))
  } else {
    stop("hypercube_vertices: dim must be greater than or equal to 1")
  }
}
