import os
import shutil
import sys

# On Windows, Python 3.8+ no longer uses PATH for DLL resolution.
# If the Fortran extension wasn't fully statically linked, we need to
# explicitly register the MinGW runtime DLL directory so that
# libgfortran, libgcc_s_seh, libquadmath, libwinpthread, etc. can be found.
if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    gfortran_path = shutil.which("gfortran")
    if gfortran_path:
        mingw_bin = os.path.dirname(os.path.realpath(gfortran_path))
        os.add_dll_directory(mingw_bin)
