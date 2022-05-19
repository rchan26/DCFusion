#' Hypercube vertices
#'
#' Obtains the vertices of the hypercube given Bessel layer information
#'
#' @param bessel_layers a list of length dim where list[[d]] is the Bessel layer
#'                      for component d, which is a list of length 4
#'                      where L denotes the hard lower bound, l denotes the soft
#'                      lower bound, u denotes the soft upper bound and U 
#'                      denotes the hard upper bound
#' @param vector vector of length dim which we check in each dimension if it
#'               occurs within the bounds - default is rep(0, dim)
#' @param transform_mat dim x dim transformation matrix which is applied to
#'                      the vertices - default is diag(1, dim)
#' @param dim dimension
#'
#' @return A list of length 2 with
#' \describe{
#'  \item{vertices}{matrix where each row is a vertex of the hypercube}
#'  \item{V}{matrix where each row is a vertex of the hypercube including points
#'           crossing the vector in each dimension}
#' }
#'
#' @export
obtain_hypercube_vertices <- function(bessel_layers,
                                      vector = rep(0, dim),
                                      transform_mat = diag(1, dim),
                                      dim) {
  if (!is.list(bessel_layers)) {
    stop("obtain_hypercube_vertices: bessel_layers must be a list of length dim")
  } else if (length(vector)!=dim) {
    stop("obtain_hypercube_vertices: vector must be a vector of length dim")
  }
  if (dim == 1) {
    if (length(bessel_layers)==1) {
      bessel_layers <- bessel_layers[[1]]
    }
    if (!identical(names(bessel_layers), c("L", "l", "u", "U"))) {
      stop("obtain_hypercube_vertices: if dim==1, bessel_layers must be a list of length 4 with names (L, l, u, U)
           or a list of length 1, where bessel_layers[[1]] is a list of length 4 with names (L, l, u, U)")
    }
    if ((bessel_layers$L < vector) & (bessel_layers$U > vector)) {
      V <- matrix(c(bessel_layers$L, vector, bessel_layers$U))
    } else {
      V <- matrix(c(bessel_layers$L, bessel_layers$U))
    }
    return(list('vertices' = matrix(c(bessel_layers$L, bessel_layers$U)) %*% transform_mat,
                'V' = V %*% transform_mat))
  } else if (dim > 1) {
    if (length(bessel_layers)!=dim) {
      stop("obtain_hypercube_vertices: if dim > 1, bessel_layers must be a list of length dim")
    } else if (!all(sapply(1:dim, function(d) identical(names(bessel_layers[[d]]), c("L", "l", "u", "U"))))) {
      stop("obtain_hypercube_vertices: if dim > 1, bessel_layers[[d]] must be a list of length 4 with names (L, l, u, U)")
    }
    bounds <- lapply(1:dim, function(d) c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
    B <- lapply(1:dim, function(d) {
      if ((bessel_layers[[d]]$L < vector[d]) & (bessel_layers[[d]]$U > vector[d])) {
        return(c(bessel_layers[[d]]$L, vector[d], bessel_layers[[d]]$U))
      } else {
        return(c(bessel_layers[[d]]$L, bessel_layers[[d]]$U))
      }
    })
    vertices <- as.matrix(expand.grid(bounds))
    V <- as.matrix(expand.grid(B))
    colnames(vertices) <- c()
    colnames(V) <- c()
    return(list('vertices' = vertices %*% transform_mat,
                'V' = V %*% transform_mat))
  } else {
    stop("obtain_hypercube_vertices: dim must be greater than or equal to 1")
  }
}
