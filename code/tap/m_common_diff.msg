1 run_model: (DF03) variable poro_unit is used before initialized
2 run_model: (DF03) variable poro_unit is used before initialized
3 run_model: (AD12) User help requested: unknown dimension 1 of array *innovuv
4 run_model: (AD12) User help requested: unknown dimension 1 of array *innovq
5 run_model: (AD12) User help requested: unknown dimension 1 of array *innovation
6 run_model: (AD12) User help requested: unknown dimension 1 of array *innovw
7 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:sum_pow_h])
8 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:qout])
9 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:zs])
10 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:av])
11 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s1])
12 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s2])
13 euler_time_step_first_b1: (AD03) Varied variable *(sporosity.phi)[_:_] written by I-O to file poro_unit
14 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
15 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s_total])
16 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:n_average])
17 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:u_mean])
18 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:v_mean])
19 calc_innovw: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:w_total])
20 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:1])])
21 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:2])])
22 update_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
23 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(fieldAccess(ident:bc, ident:sum_mass_flux), expressions[ident:num_bc])])
24 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(pointerAccess(fieldAccess(ident:bc, ident:sum_mass_flux), none()), expressions[ident:num_bc])])
25 default_values: (TC16) Type mismatch in assignment: INTEGER*4 receives REAL*4
26 alloc_or_realloc_station: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
27 alloc_or_realloc_station: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
28 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
29 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
30 end_mpi: (TC32) Conflicting numbers of arguments for procedure exit, expected 1, is here 0
31 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
32 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
33 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
34 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
35 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
36 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
37 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
38 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
39 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
40 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
41 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
42 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
43 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
44 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
45 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
46 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
47 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
48 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
49 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
50 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
51 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
52 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
53 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
54 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
55 reading_args: (TC33) External procedure command_argument_count is not declared
56 reading_args: (TC33) External procedure get_command_argument is not declared
57 msh_initialise: (TC33) External procedure mesh_input is not declared
58 mpi_sum_r: Tool: (Build MP call info) Unexpected MPI primitive mpi_wait_all in call(none(), ident:mpi_wait_all, expressions[])
