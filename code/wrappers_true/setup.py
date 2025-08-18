from setuptools import setup

setup(
    name="dassflow2d",
    version="1.0",    description="SWE modeling, variational data assimilation, parallel calculus",
    url="/",
    author="INRAE/INSA",
    
    packages=["dassflow2d", 
    	       "dassflow2d.wrapping", 
    	       "dassflow2d.core", 
    	       "dassflow2d.assim"],
    	       
    install_requires=[
        "matplotlib>=3.4.1", # classical plot library
        "numpy>=1.21.2",    # array mannipulation
        "pandas>=1.3.5",    # data manipulation,  not used yet
        "scipy>=1.7.1" ,   # minimization
        "h5py>=3.6.0"  ,  # save as hdf5
        "pyvista"      ,  # mesh visu and manipulation in python
        "vtk"          ,   # pyvista use vtk
       # "f90wrap"         # for wrapping
        # autodoc
    ],
    package_data={
        "dassflow2d": ["wrapping/_wrapping*.so"]
    },
    include_package_data=True,
    zip_safe=False,
)
