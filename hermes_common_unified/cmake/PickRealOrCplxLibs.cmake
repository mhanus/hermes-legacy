#
# PETSc and MUMPS have different set of libraries for real and complex versions.
# PETSc also has specific include directories for real and complex versions, 
# apart from the architecture-independent ones. Following macros decide which 
# versions to use, according to which version of $HERMES we link the current 
# target $TRGT against. 
#
# Note that we cannot use 'include_directories' to add the architecture-specific
# include directories for PETSc as this function does not distinguish between
# build targets - these directories are thus excluded from dependency scanning, 
# but this is not an issue as long as there aren't any other libraries with
# same-named include files.
#
# Unfortunately, the latter happens with the dummy sequential "mpi.h" file, which
# is used (with this name) both in PETSc and MUMPS. Therefore, it is neccessary 
# to treat this single file via include_directories.
#
macro(SET_PETSC_FLAGS TRGT ARSP_INCLUDE_DIRS)
  if(WITH_PETSC)
    get_property(ARCH_SPECIFIC_FLAGS TARGET ${TRGT} PROPERTY COMPILE_FLAGS)
    foreach(_DIR ${${ARSP_INCLUDE_DIRS}})
      set(ARCH_SPECIFIC_FLAGS "${ARCH_SPECIFIC_FLAGS} -I${_DIR}")
    endforeach(_DIR)
    set_property(TARGET ${TRGT} PROPERTY COMPILE_FLAGS ${ARCH_SPECIFIC_FLAGS})
  endif(WITH_PETSC)
endmacro(SET_PETSC_FLAGS)

macro(PICK_REAL_OR_CPLX_LIBS HERMES_COMMON TRGT)
    set(PETSC_LIBRARIES ${PETSC_CPLX_LIBRARIES})
    set(MUMPS_LIBRARIES ${MUMPS_CPLX_LIBRARIES})
    LIST(APPEND MUMPS_LIBRARIES ${MUMPS_REAL_LIBRARIES})
endmacro(PICK_REAL_OR_CPLX_LIBS)    

macro(PICK_REAL_OR_CPLX_INCS HERMES_COMMON TRGT)
#we use only complex petsc (even for computing in real numbers)
    SET_PETSC_FLAGS(${TRGT} PETSC_CPLX_INCLUDE_DIRS)
endmacro(PICK_REAL_OR_CPLX_INCS HERMES_COMMON TRGT)  
