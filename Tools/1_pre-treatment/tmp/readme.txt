f90wrap -m gen_case _gen_channel_case.f90 --package
f2py-f90wrap -c -m gen_case _gen_channel_case.f90 f90wrap_toplevel.f90
# compile wraped code to generate input files
python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90
# generate input file
python3 run.py


# to  macdonald case
python3 gen_macdo.py
