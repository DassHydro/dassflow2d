###############################################
###############################################
# RENAME VTK FILE SO THAT THEY CAN BE USED AS A "TIME SERIE" IN PARAVIEW
###############################################
###############################################
source_dir <- "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/res"
  #"/home/lilian.villenave/Documents/save/GIT/V2_dassflow/dassflow2D-devlilian/trunk/code/bin_AUDE_unif_30m_init/res"
write_dir <- file.path(source_dir, "vtk_2")

dir.create(write_dir)
name_files <- list.files(source_dir)
path_files <- list.files(source_dir, full.names = TRUE)

vtk_id <- grepl(pattern = ".vtk", x = name_files)

path_files <- path_files[vtk_id]
name_files <- name_files[vtk_id]

time <- sapply(strsplit(name_files, "_"), function(x) return(x[2]))
time <- gsub(x = time, pattern = ".vtk", replacement = "")
time <- gsub(x = time, pattern = "initial", replacement = "0")
time <- gsub(x = time, pattern = "final", replacement = "1E+010")
new_time <- as.numeric(time)
new_id <- order(new_time)

path_files <- path_files[new_id]
new_path_files <-  paste0("dassflow_res_", stringr::str_pad(seq(0,length(path_files)-1), 5, pad = "0"), ".vtk")


file.copy(from = path_files, to = file.path(write_dir,new_path_files), overwrite = TRUE)

write.csv(x= cbind(path_files, new_path_files), file = file.path(write_dir, "correspondance_timestep.csv"))



