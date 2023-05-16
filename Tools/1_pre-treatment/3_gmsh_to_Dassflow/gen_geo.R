tmp = readLines("/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/tmp9.msh")
new_file = "/home/livillenave/Documents/software/gmsh-4.10.5-Linux64/lilian.geo"


all_info_node = tmp[ grep("Nodes",tmp)[1]+1]
nb_node = as.numeric(strsplit(all_info_node, " ")[[1]][2])
start_node = grep("Nodes",tmp)[1]+2
end_node = grep("Nodes",tmp)[2]-1

id_cells = grep("Elements",tmp)
elements_metadata =tmp[id_cells+1]
elements_metadata = strsplit(elements_metadata[1], " ")[[1]]
unkown_elements = as.integer(elements_metadata[1])
nb_cell = as.integer(elements_metadata[2])
start_c = 1974#id_cells[1]+1+unkown_elements
end_c = id_cells[2]
nb_cell_lines =  end_c -start_c 

cat( "# Generated mesh with R script gen_geo.R ||| number of nodes |   number of cells | mesh scale == 0 always",
    file = new_file,
    append = FALSE)

cat( paste0("\n",nb_node," ",nb_cell_lines , "  0 \n"),
     file = new_file,
     append = TRUE)

cat( "#Nodes||| id node, x coord, y coord, bathy (x,y) \n",
     file = new_file,
     append = TRUE)



all_nodes_data = strsplit(tmp[start_node:end_node], " ")
len <- sapply(all_nodes_data, function(x){return(length(x))} )
lines_nodedata= which(len == 3)
id_node = 0


all_cell_data = strsplit(tmp[start_c:end_c], " ")


for( i in lines_nodedata){
  
  id_node = id_node + 1
  
  cat( 
  paste(c(format(id_node, scientific=F), all_nodes_data[[i]],"\n"), collapse = " "),
  file = new_file,
  append = TRUE  )
  
}


cat( "# Cells|||  id of cell, id node 1, id node 2, id  node 3,  id node 4 , land_type, bathymetry \n",
     file = new_file,
     append = TRUE)

for( i in 1:nb_cell_lines){
  to_paste = all_cell_data[[i]]
  res = c(i, to_paste[c(2,3,4)], to_paste[2], 1, 0,"\n")
  cat( res,
    file = new_file,
    append = TRUE  )
}



cat( "# Boundaries \n",
     file = new_file,
     append = TRUE)

cat( "INLET 0 0 \n",
     file = new_file,
     append = TRUE)

cat( "OUTLET 0 0 \n",
     file = new_file,
     append = TRUE)





