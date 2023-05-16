from setuptools import setup

setup(
    name="gen_case",
    version="0.1",
    description="test",
    url="/",
    author="INRAE",
    
    packages=["gen_case"],
    	       
    install_requires=[
        "matplotlib>=3.4.1", # classical plot library
        "numpy>=1.21.2",    # array mannipulation
        "pandas>=1.3.5",    # data manipulation,  not used yet
        "scipy>=1.7.1" ,   # minimization
        "h5py>=3.6.0"  ,  # save as hdf5
        "pyvista"      ,  # mesh visu and manipulation in python
        "vtk"          ,   # pyvista use vtkS
    ],
    package_data={
        "gen_case": ["gen_case*.so"]
    },
    include_package_data=True,
    zip_safe=False,
)
