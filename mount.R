mount_homothetic <- function(a,tDeg) {
  bind_rows(get_ellipse_point(a,tDeg),
            get_ellipse_point(a,tDeg+120),
            get_ellipse_point(a,tDeg-120))
}

mount_fs <- function(a,tDeg) {
  c <- sqrt(a*a-1)
  get_ellipse_foci(a) %>%
    bind_rows(get_ellipse_point(a,tDeg))
}

mount_major <- function(a,tDeg) {
  tibble(x=c(-a,a),y=c(0,0)) %>%
    bind_rows(get_ellipse_point(a,tDeg))
}
