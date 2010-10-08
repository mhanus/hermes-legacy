#
# PETSc
#
# Only a single configuration is currently supported in FEMHUB (the real 
# sequential one by default). It is also limited to PETSc 3.1, which is the first
# version that creates one file for the whole library by default.
#
SET(MY_PETSC_LIB_DIRS $ENV{MY_PETSC_LIB_DIRS})
SET(MY_PETSC_INC_DIRS $ENV{MY_PETSC_INC_DIRS})

IF(NOT MY_PETSC_LIB_DIRS OR NOT MY_PETSC_INC_DIRS)
  SET(MY_PETSC_LIB_DIRS $ENV{SAGE_LOCAL}/lib)
  SET(MY_PETSC_INC_DIRS $ENV{SAGE_LOCAL}/include)
ENDIF(NOT MY_PETSC_LIB_DIRS OR NOT MY_PETSC_INC_DIRS)

# Try to find petsc.h in the root include directory.
FIND_PATH(COMMON_PETSC_INCLUDE_DIRS petsc.h PATHS ${MY_PETSC_INC_DIRS} NO_DEFAULT_PATH)

IF(COMMON_PETSC_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${MY_PETSC_INC_DIRS}/mpiuni)
  # Add path for the sequential emulation of MPI.
  SET(COMMON_PETSC_INCLUDE_DIRS ${COMMON_PETSC_INCLUDE_DIRS} ${MY_PETSC_INC_DIRS}/mpiuni)
ENDIF(COMMON_PETSC_INCLUDE_DIRS AND NOT WITH_MPI AND EXISTS ${MY_PETSC_INC_DIRS}/mpiuni)

# In Femhub, there is either the real, or the complex lib (not both together).    

FIND_LIBRARY(PETSC_LIB petsc ${MY_PETSC_LIB_DIRS}  NO_DEFAULT_PATH)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG PETSC_LIB COMMON_PETSC_INCLUDE_DIRS)
   
# linux specific (?)
SET(${PETSC_LIB} ${PETSC_LIB} dl)
    
SET(PETSC_REAL_LIBRARIES ${PETSC_LIB})     
SET(PETSC_CPLX_LIBRARIES ${PETSC_LIB})          
