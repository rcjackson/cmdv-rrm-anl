Analysis of 3D winds with Python code
=====================================

In this directory are codes that are used to process multidop grids in parallel as shown on Robert Jackson's SciPy 2017 Poster.  
Dependencies:
    multidop (available at: https://github.com/nasa/MultiDop) 
    
    Python ARM Radar Toolkit (available at: https://github.com/ARM-DOE/pyart)
    
    parse (available at: https://github.com/r1chardj0n3s/parse)
    
    netcdf4
    
    xarray
    
    dask
    
    distributed
    
Brief description of each file:

**time_procedures.py** - This module contains code that catalogs a list of radar and sounding files assuming that the timestamp is in the name of the file. In order to use this for your radar dataset, the file_name_str and format_str files need to be modified to match your dataset's naming convention.

**multidop_run.py** - This module contains the code that calls dask to run multidop on multiple scans in parallel. 

**multidop_time.py** - This module contains the code that retrieved the 3D wind grids from two radars with dealiased velocity data.

**w_pdfs_by_drosdowsky.py** - This code provides an example of how dask is used to derive p.d.f.s of vertical velocity from thousands of Cf-Complaint grids that are produced by multidop. 
                     
