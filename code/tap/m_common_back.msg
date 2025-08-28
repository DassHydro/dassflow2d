1 run_model: (DF03) variable poro_unit is used before initialized
2 run_model: (DF03) variable poro_unit is used before initialized
3 run_model: (AD12) User help requested: unknown dimension 1 of array *dof%infil
4 run_model: (AD12) User help requested: unknown dimension 1 of array *dof%v
5 run_model: (AD12) User help requested: unknown dimension 1 of array *dof%u
6 run_model: (AD12) User help requested: unknown dimension 1 of array *dof%h
7 run_model: (AD12) User help requested: unknown dimension 1 of array *innovuv
8 run_model: (AD12) User help requested: unknown dimension 1 of array *innovq
9 run_model: (AD12) User help requested: unknown dimension 1 of array *innovation
10 run_model: (AD12) User help requested: unknown dimension 1 of array *sporosity%a
11 run_model: (AD12) User help requested: unknown dimension 1 of array *sporosity%phi
12 run_model: (AD12) User help requested: unknown dimension 1 of array *bathy_cell
13 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%sum_mass_flux
14 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
15 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
16 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rat
17 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%outflow
18 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%inflow
19 run_model: (AD12) User help requested: unknown dimension 1 of array *ptf
20 run_model: (AD12) User help requested: unknown dimension 1 of array *ptf
21 run_model: (AD12) User help requested: unknown dimension 1 of array *ptf
22 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
23 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
24 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rat
25 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rat
26 run_model: (AD12) User help requested: unknown dimension 1 of array *bc%hyd
27 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
28 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rat
29 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%outflow
30 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%inflow
31 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%h
32 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%h
33 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%sum_mass_flux
34 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *sporosity%a
35 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *sporosity%phi
36 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bathy_cell
37 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%sum_mass_flux
38 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%rain
39 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%h
40 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%h
41 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%u
42 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%u
43 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%v
44 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%v
45 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%infil
46 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *dof%infil
47 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *bc%sum_mass_flux
48 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *innovation
49 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *innovuv
50 sub_run_model: (AD12) User help requested: unknown dimension 1 of array *innovq
51 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:sum_pow_h])
52 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:qout])
53 set_bc: (AD07) Data-Flow recovery (TBR) needs to restore allocatable pointer *(bc.inflow), not implemented yet
54 set_bc: (AD12) User help requested: unknown dimension 1 of array *bc%inflow
55 set_bc: (AD07) Data-Flow recovery (TBR) needs to restore allocatable pointer *(bc.outflow), not implemented yet
56 set_bc: (AD12) User help requested: unknown dimension 1 of array *bc%outflow
57 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:zs])
58 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:av])
59 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s1])
60 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s2])
61 euler_time_step_first_b1: (AD03) Varied variable *(sporosity.phi)[_:_] written by I-O to file poro_unit
62 euler_time_step_first_b1: (AD12) User help requested: unknown dimension 1 of array *bathy_cell
63 euler_time_step_first_b1: (AD12) User help requested: unknown dimension 1 of array *bc%sum_mass_flux
64 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
65 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s_total])
66 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:n_average])
67 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:u_mean])
68 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:v_mean])
69 calc_innovw: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:w_total])
70 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:1])])
71 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:2])])
72 calc_cost_function: (AD12) User help requested: unknown dimension 1 of array *bc%hyd
73 calc_cost_function: (AD12) User help requested: unknown dimension 1 of array *bc%hyd
74 update_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
75 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(fieldAccess(ident:bc, ident:sum_mass_flux), expressions[ident:num_bc])])
76 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(pointerAccess(fieldAccess(ident:bc, ident:sum_mass_flux), none()), expressions[ident:num_bc])])
77 default_values: (TC16) Type mismatch in assignment: INTEGER*4 receives REAL*4
78 alloc_or_realloc_station: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
79 alloc_or_realloc_station: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
80 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
81 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
82 end_mpi: (TC32) Conflicting numbers of arguments for procedure exit, expected 1, is here 0
83 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
84 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
85 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
86 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
87 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
88 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
89 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
90 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
91 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
92 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
93 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
94 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
95 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
96 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
97 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
98 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
99 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
100 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
101 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
102 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
103 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
104 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
105 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
106 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
107 reading_args: (TC33) External procedure command_argument_count is not declared
108 reading_args: (TC33) External procedure get_command_argument is not declared
109 msh_initialise: (TC33) External procedure mesh_input is not declared
110 mpi_sum_r: Tool: (Build MP call info) Unexpected MPI primitive mpi_wait_all in call(none(), ident:mpi_wait_all, expressions[])
