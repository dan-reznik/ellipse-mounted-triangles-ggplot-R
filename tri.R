law_of_cosines <- function(a,b,c) { (b*b+c*c-a*a)/(2*b*c) }

tri_dist <- function(df_tri,v1,v2) {
  edist(df_tri[[v1,1]], df_tri[[v1,2]],
        df_tri[[v2,1]], df_tri[[v2,2]])
}

tri_sides <- function(df_tri) {
  s1 <- tri_dist(df_tri,2,3)
  s2 <- tri_dist(df_tri,3,1)
  s3 <- tri_dist(df_tri,1,2)
  c(s1,s2,s3)
}

tri_cosines <- function(sides) {
  a <- sides[1]; b <- sides[2]; c <- sides[3];
  x <- law_of_cosines(a,b,c)
  y <- law_of_cosines(b,c,a)
  z <- law_of_cosines(c,a,b)
  c(x,y,z)
}

dot3 <- function(v1,v2) {
  v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]
}

trilin_to_cartesian <- function(df_tri,sides,ts) {
  v <- map2_dbl(sides,ts,~.x*.y)
  # c(sides[1]*ts[1],sides[2]*ts[2],sides[3]*ts[3])
  denom <- sum(v) # v[1]+v[2]+v[3]
  x <- dot3(v, df_tri$x) / denom
  y <- dot3(v, df_tri$y) / denom 
  tibble(x=x,y=y)
}

tri_generic <- function(df_tri, sides, ts_list) {
  ts_list %>%
    map_dfr(~trilin_to_cartesian(df_tri,sides,.x))
}

tri_orthic <- function(df_tri) {
  sides <- tri_sides(df_tri)
  cs <- tri_cosines(sides)
  secs <- cs %>% map_dbl(~1/.x)
  ts_list <- list(row1=c(0,secs[2],secs[3]),
                  row2=c(secs[1],0,secs[3]),
                  row3=c(secs[1],secs[2],0)) 
  tri_generic(df_tri,sides,ts_list)
}

# https://faculty.evansville.edu/ck6/encyclopedia/ETC.html
# https://github.com/dan-reznik/ellipse-mounted-loci-p5js
tri_X131 <- function(df_tri) {
  sides <- tri_sides(df_tri)
  a <- sides[1]; b <- sides[2]; c <- sides[3]
  c2 <- c*c; c4 <- c2*c2; c6 <- c2*c4
  b2 <- b*b; b4 <- b2*b2; b6 <- b2*b4
  a2 <- a*a; a4 <- a2*a2; a6 <- a2*a4
  c8 <- c2*c6; b8 <- b2*b6; a8 <- a2*a6
  v1 <- (a2-b2-c2)*(a4*b2-2*a2*b4+b6+a4*c2+2*a2*b2*c2-b4*c2-2*a2*c4-b2*c4+c6)*(2*a8-3*a6*b2+a4*b4-a2*b6+b8-3*a6*c2+2*a4*b2*c2+a2*b4*c2-4*b6*c2+a4*c4+a2*b2*c4+6*b4*c4-a2*c6-4*b2*c6+c8)
  v2 <- (-a2+b2-c2)*(a6-2*a4*b2+a2*b4-a4*c2+2*a2*b2*c2+b4*c2-a2*c4-2*b2*c4+c6)*(a8-a6*b2+a4*b4-3*a2*b6+2*b8-4*a6*c2+a4*b2*c2+2*a2*b4*c2-3*b6*c2+6*a4*c4+a2*b2*c4+b4*c4-4*a2*c6-b2*c6+c8)
  v3 <- (-a2-b2+c2)*(a6-a4*b2-a2*b4+b6-2*a4*c2+2*a2*b2*c2-2*b4*c2+a2*c4+b2*c4)*(a8-4*a6*b2+6*a4*b4-4*a2*b6+b8-a6*c2+a4*b2*c2+a2*b4*c2-b6*c2+a4*c4+2*a2*b2*c4+b4*c4-3*a2*c6-3*b2*c6+2*c8)
  ts <- c(v1/a,v2/b,v3/c)
  trilin_to_cartesian(df_tri, sides, ts)
}

tri_X138 <- function(df_tri) {
  sides <- tri_sides(df_tri)
  a <- sides[1]; b <- sides[2]; c <- sides[3]
  c2 <- c*c; c4 <- c2*c2; c6 <- c2*c4
  b2 <- b*b; b4 <- b2*b2; b6 <- b2*b4
  a2 <- a*a; a4 <- a2*a2; a6 <- a2*a4
  c8 <- c2*c6; b8 <- b2*b6; a8 <- a2*a6
  v1 <- (a2+b2-c2)*(a2-b2+c2)*(a2*b2-b4+a2*c2+2*b2*c2-c4)*(2*a8-4*a6*b2+a4*b4+2*a2*b6-b8-4*a6*c2+4*a4*b2*c2-2*a2*b4*c2+2*b6*c2+a4*c4-2*a2*b2*c4-2*b4*c4+2*a2*c6+2*b2*c6-c8)*(a8-2*a6*b2+2*a4*b4-2*a2*b6+b8-2*a6*c2-a4*b2*c2+2*a2*b4*c2+b6*c2+2*a4*c4+2*a2*b2*c4-4*b4*c4-2*a2*c6+b2*c6+c8);
  v2 <- (a2+b2-c2)*(-a2+b2+c2)*(-a4+a2*b2+2*a2*c2+b2*c2-c4)*(-a8+2*a6*b2+a4*b4-4*a2*b6+2*b8+2*a6*c2-2*a4*b2*c2+4*a2*b4*c2-4*b6*c2-2*a4*c4-2*a2*b2*c4+b4*c4+2*a2*c6+2*b2*c6-c8)*(a8-2*a6*b2+2*a4*b4-2*a2*b6+b8+a6*c2+2*a4*b2*c2-a2*b4*c2-2*b6*c2-4*a4*c4+2*a2*b2*c4+2*b4*c4+a2*c6-2*b2*c6+c8);
  v3 <- (a2-b2+c2)*(-a2+b2+c2)*(-a4+2*a2*b2-b4+a2*c2+b2*c2)*(a8+a6*b2-4*a4*b4+a2*b6+b8-2*a6*c2+2*a4*b2*c2+2*a2*b4*c2-2*b6*c2+2*a4*c4-a2*b2*c4+2*b4*c4-2*a2*c6-2*b2*c6+c8)*(-a8+2*a6*b2-2*a4*b4+2*a2*b6-b8+2*a6*c2-2*a4*b2*c2-2*a2*b4*c2+2*b6*c2+a4*c4+4*a2*b2*c4+b4*c4-4*a2*c6-4*b2*c6+2*c8);
  ts <- c(v1/a,v2/b,v3/c)
  trilin_to_cartesian(df_tri, sides, ts)
}