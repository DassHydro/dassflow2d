library(raster)
library(rgdal)

name_output = "aude_v1.geo"


#shape_majeur = readOGR( "/home/livillenave/Images/test_again.shp")
shape_majeur = readOGR( "/home/livillenave/Documents/distant/SD-FLOOD/real_case-AUDE/DATA/DASSFLOW/V2/hydraulic/tmp/contour1.shp")
points_majeur =  shape_majeur@polygons[[1]]@Polygons[[1]]@coords #shape_majeur@lines[[1]]@Lines[[1]]@coords
nbpoint_majeur = nrow(points_majeur)
plot(shape_majeur)

#shape_mineur = readOGR("/home/livillenave/lit_mineur_crop.shp")
shape_mineur = readOGR("/home/livillenave/Documents/distant/SD-FLOOD/real_case-AUDE/DATA/DASSFLOW/V2/hydraulic/tmp/contour2.shp")
points_mineur = shape_mineur@polygons[[1]]@Polygons[[1]]@coords # shape_mineur@lines[[1]]@Lines[[1]]@coords
nbpoint_mineur = nrow(points_mineur)
plot(shape_mineur, add = TRUE, col = "green")



cat("lc_majeur=200;",
    file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output) ,
    append = FALSE)

cat("\nlc_mineur=200;",
    file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output) ,
    append = TRUE)

#========================================================================#
# GENERATE points
#========================================================================#

# ------------------ LIT MAJEUR --------------- #
for(i in 1:nbpoint_majeur){
  cat(paste0("\nPoint(",i, ") = {", points_majeur[i,1], ",", points_majeur[i,2], ", 0 ,", "lc_majeur};"),
      file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/", name_output),
      append = TRUE)
}

# ------------------ LIT MINEUR --------------- #
for(i in 1:nbpoint_mineur){
  cat(paste0("\nPoint(",i+nbpoint_majeur, ") = {", points_mineur[i,1], ",", points_mineur[i,2], ", 0 ,", "lc_mineur};"),
      file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64",name_output),
      append = TRUE)
}

#========================================================================#
# GENERATE lines
#========================================================================#

# ------------------ LIT MAJEUR --------------- #
for(i in 1:(nbpoint_majeur-1)){
  if(i !=nbpoint_majeur-1){
  cat(paste0("\nLine(",i, ") = {", i, ",", i+1, "};"),
      file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
      append = TRUE)
  }else{
    cat(paste0("\nLine(",i, ") = {", i, ",", 1, "};"),
    file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
    append = TRUE)
  }
}

# ------------------ LIT MINEUR --------------- #

for(i in 1:(nbpoint_mineur-1)){
  if(i !=nbpoint_mineur-1){
    cat(paste0("\nLine(",nbpoint_majeur+i, ") = {", nbpoint_majeur+i, ",", nbpoint_majeur+i+1, "};"),
        file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
        append = TRUE)
  }else{
    cat(paste0("\nLine(",nbpoint_majeur+i, ") = {",nbpoint_majeur+ i, ",",nbpoint_majeur+ 1, "};"),
        file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
        append = TRUE)
  }
}

#========================================================================#
# CURVE LOOP + meshing surface
#========================================================================#

cat(paste0("\nCurve Loop(1) = {1:",nbpoint_majeur-1, "};"),
    file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
    append = TRUE)

cat("\nPlane Surface(1) = {1};",
    file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
    append = TRUE)


# ! add lit_mineur specification
for(i in 1:(nbpoint_mineur-1)){
  cat(paste0("\nLine{",nbpoint_majeur+i, "} In Surface{1};"),
      file = file.path("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/",name_output),
      append = TRUE)
}

