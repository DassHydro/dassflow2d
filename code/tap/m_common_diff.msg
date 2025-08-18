1 run_model: (AD12) User help requested: unknown dimension 1 of array *mesh%edge
2 run_model: (AD12) User help requested: unknown dimension 1 of array *mesh%edge
3 run_model: (AD12) User help requested: unknown dimension 1 of array *innovuv
4 run_model: (AD12) User help requested: unknown dimension 1 of array *innovation
5 run_model: (AD12) User help requested: unknown dimension 1 of array *innovw
6 run_model: (AD12) User help requested: unknown dimension 1 of array *innovq
7 sub_run_model: (TC33) External procedure euler_time_step_first_b1 is not declared
8 sub_run_model: (DF05) External routine euler_time_step_first_b1, called by (sub_run_model), is assumed default In:..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9)..........(0)..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9)..........(0)..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9).111111111(0)1111111111(1)1111111111(2)1111111111(3)1111111111(4)1111111111(5)1111111111(6)1111111111(7)1111111111(8)1 and Out:..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9)..........(0)..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9)..........(0)..........(1)..........(2)..........(3)..........(4)..........(5)..........(6)..........(7)..........(8)..........(9).111111111(0)1111111111(1)1111111111(2)1111111111(3)1111111111(4)1111111111(5)1111111111(6)1111111111(7)1111111111(8)1
9 sub_run_model: (AD09) Please provide a differential of function euler_time_step_first_b1 for arguments Arg1:*(_.h)[_:_]=(in;out) Arg1:*(_.u)[_:_]=(in;out) Arg1:*(_.v)[_:_]=(in;out) Arg1:*(_.infil)[_:_]=(in;out) Arg1:_.t_display=(in;out) Arg1:*(_.grad_h)[_:_].x=(in;out) Arg1:*(_.grad_h)[_:_].y=(in;out) Arg1:*(_.grad_u)[_:_].x=(in;out) Arg1:*(_.grad_u)[_:_].y=(in;out) Arg1:*(_.grad_v)[_:_].x=(in;out) Arg1:*(_.grad_v)[_:_].y=(in;out) Arg1:*(_.grad_z)[_:_].x=(in;out) Arg1:*(_.grad_z)[_:_].y=(in;out) Arg2:*(_.node)[_:_].coord.x=(in;out) Arg2:*(_.node)[_:_].coord.y=(in;out) Arg2:*(_.cell)[_:_].surf=(in;out) Arg2:*(_.cell)[_:_].invsurf=(in;out) Arg2:*(_.cell)[_:_].peri=(in;out) Arg2:*(_.cell)[_:_].grav.x=(in;out) Arg2:*(_.cell)[_:_].grav.y=(in;out) Arg2:*(_.cellb)[_:_].grav.x=(in;out) Arg2:*(_.cellb)[_:_].grav.y=(in;out) Arg2:*(_.edge)[_:_].length=(in;out) Arg2:*(_.edge)[_:_].center.x=(in;out) Arg2:*(_.edge)[_:_].center.y=(in;out) Arg2:*(_.edge)[_:_].normal.x=(in;out) Arg2:*(_.edge)[_:_].normal.y=(in;out) Arg2:*(_.edge)[_:_].tangent.x=(in;out) Arg2:*(_.edge)[_:_].tangent.y=(in;out) Arg2:*(_.edge)[_:_].vcell.x=(in;out) Arg2:*(_.edge)[_:_].vcell.y=(in;out) Arg2:*(_.edge)[_:_].v_edge_cell[1:2].x=(in;out) Arg2:*(_.edge)[_:_].v_edge_cell[1:2].y=(in;out) Arg2:_.scal=(in;out) Arg2:_.surf=(in;out)
10 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:sum_pow_h])
11 set_bc: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:qout])
12 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:zs])
13 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:av])
14 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s1])
15 newton_qin: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s2])
16 euler_time_step_first_b1_poro: (TC33) External procedure update_all_porosities is not declared
17 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
18 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:s_total])
19 calc_innovation: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_i in call(none(), ident:mpi_sum_i, expressions[ident:n_average])
20 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:u_mean])
21 calc_innovuv: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:v_mean])
22 calc_innovw: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:w_total])
23 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:1])])
24 calc_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(ident:cost_part, expressions[intCst:2])])
25 update_cost_function: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[ident:h_mean])
26 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(fieldAccess(ident:bc, ident:sum_mass_flux), expressions[ident:num_bc])])
27 sw_post_treatment: Tool: (Build MP call info) Unexpected MPI primitive mpi_sum_r in call(none(), ident:mpi_sum_r, expressions[arrayAccess(pointerAccess(fieldAccess(ident:bc, ident:sum_mass_flux), none()), expressions[ident:num_bc])])
28 default_values: (TC16) Type mismatch in assignment: INTEGER*4 receives REAL*4
29 alloc_or_realloc_station: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
30 alloc_or_realloc_station: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBS(*) POINTER
31 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
32 alloc_or_realloc_stationq: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE STATION_OBSQ(*) POINTER
33 end_mpi: (TC32) Conflicting numbers of arguments for procedure exit, expected 1, is here 0
34 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
35 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
36 reallocate_cell: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
37 reallocate_cell: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPE(*) POINTER
38 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
39 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
40 reallocate_edge: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
41 reallocate_edge: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPE(*) POINTER
42 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
43 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
44 reallocate_node: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
45 reallocate_node: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPE(*) POINTER
46 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
47 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
48 reallocate_cellb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
49 reallocate_cellb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE CELLTYPELIM(*) POINTER
50 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
51 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
52 reallocate_edgeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
53 reallocate_edgeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE EDGETYPELIM(*) POINTER
54 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
55 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
56 reallocate_nodeb: (TC30) Type mismatch in argument 1 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
57 reallocate_nodeb: (TC30) Type mismatch in argument 2 of procedure move_alloc, expected REAL*8(*) POINTER, is here derived TYPE NODETYPELIM(*) POINTER
58 reading_args: (TC33) External procedure command_argument_count is not declared
59 reading_args: (TC33) External procedure get_command_argument is not declared
60 msh_initialise: (TC33) External procedure mesh_input is not declared
61 mpi_sum_r: Tool: (Build MP call info) Unexpected MPI primitive mpi_wait_all in call(none(), ident:mpi_wait_all, expressions[])
