toRad <- function(tDeg) { tDeg*pi/180.0 }
edist <- function(x1,y1,x2,y2) {
  dx <- x1-x2
  dy <- y1-y2
  sqrt(dx*dx+dy*dy)
}

get_ellipse_point <- function(a,tDeg) {
  t <- toRad(tDeg)
  tibble(x=a*sin(t),y=cos(t))
}

get_ellipse_foci <- function(a) {
  c = sqrt(a*a-1)
  tibble(x=c(-c,c),y=c(0,0))
}

get_locus <- function(a,degStep,tri_fn,tri_derived_fn,x_fn) {
  degs <- seq(0,360,degStep)
  tris <- degs %>% map(~{tri_fn(a,.x)})
  tris_derived <- tris %>% map(tri_derived_fn)
  df_xs <- tris_derived %>% map_dfr(x_fn)
  df_xs %>% mutate(deg=degs,.before=x)
}

draw_scene <- function(a,tDeg,tri_fn,tri_derived_fn,X_fn,
                       df_locus) {
  df_tri <- tri_fn(a,tDeg)
  df_derived <- tri_derived_fn(df_tri) 
  df_foci <- get_ellipse_foci(a)
  df_X <- X_fn(df_derived)
  ggplot() +
    geom_ellipse(aes(a = a, b = 1,x0=0,y0=0,angle=0), size=1.0,color="black") +
    geom_point(aes(x,y),color="black",data=df_foci) +
    geom_polygon(aes(x,y),fill=NA, size=1.25,color="blue", data=df_tri) + 
    geom_point(aes(x,y),color="blue",size=2.5, data=df_derived) +
    geom_polygon(aes(x,y),fill=NA, size=1.25,color="orange", data=df_derived) + 
    geom_point(aes(x,y),color="orange",size=2.5, data=df_derived) +
    geom_path(aes(x,y),color="red",size=1.25,data=df_locus) +
    geom_point(aes(x,y),color="red",size=2.5, data=df_X) +
    geom_text(aes(x,y),label="X131",color="red",data=df_X,vjust=-.2,hjust=-.2) +
    coord_fixed() +
    theme_minimal()
}

calc_and_draw <- function(a,tDeg,tDegStep,
                          tri_fn,tri_derived_fn,X_fn) {
  df_locus <- get_locus(a,tDegStep,tri_fn,
                        tri_derived_fn,X_fn)
  draw_scene(a,tDeg,tri_fn,tri_derived_fn,X_fn,df_locus)
}
