f90wrap -m module _gen_channel_case.f90 --package
# compile wraped code to generate input files
python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90
# generate input file
python3 run.py


# to  macdonald case
python3 gen_macdo.py
