#
# PARDISO
#
# set WITH_PARDISO to YES to enable pardiso support
# set WITH_OPENMP to YES if using the parallel version of the library
#
# PARDISO requires the user to agree with a special license agreement, manually download a 
# precompiled library for his architecture, let a license number be generated for him and finally
# store this number in file 'pardiso.lic' (required location of which is described in chapter 3 of
# PARDISO manual). There is no way to automate this process in femhub, so in order to use PARDISO,
# the user must ensure that all requirements are met and point to the library file in the 
# environment variable MY_PARDISO_LIB before compiling femhub (e.g. by saying 
# 'export MY_PARDISO_LIB=/opt/pardiso/libpardiso400_GNU430_AMD_IA64.so').
#

IF(NOT EXISTS $ENV{MY_PARDISO_LIB})
  SET(PARDISO_FOUND FALSE)	
ELSE(NOT EXISTS $ENV{MY_PARDISO_LIB})
  SET(PARDISO_LIBRARY $ENV{MY_PARDISO_LIB})
  SET(PARDISO_FOUND TRUE)
ENDIF(NOT EXISTS $ENV{MY_PARDISO_LIB})

IF (PARDISO_FOUND)
	IF (NOT PARDISO_FIND_QUIETLY)
		MESSAGE(STATUS "Found PARDISO: ${PARDISO_LIBRARY}")
	ENDIF (NOT PARDISO_FIND_QUIETLY)
ELSE (PARDISO_FOUND)
	IF (PARDISO_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find PARDISO. You must specify full path to the library in the environment variable MY_PARDISO_LIB, e.g. 'export MY_PARDISO_LIB=/opt/pardiso/libpardiso400_GNU430_AMD_IA64.so'.")
	ENDIF (PARDISO_FIND_REQUIRED)
ENDIF (PARDISO_FOUND)
