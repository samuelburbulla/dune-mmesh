import os
path = os.path.join( os.path.dirname(__file__), "tutorial" )
execute  = "cp -r " + path + " "
execute += "mmesh_tutorial"
status = os.system(execute)
if status != 0: raise RuntimeError(status)

print("##########################################################################")
print("## Some example scripts are now located in the 'mmesh_tutorial' folder. ##")
try:
    import gmsh
except ImportError:
    print("## Note: the examples requires the installation of 'gmsh'.")
try:
    import matplotlib
except ImportError:
    print("## Note: the examples requires the installation of 'matplotlib'.")
print("##########################################################################")
