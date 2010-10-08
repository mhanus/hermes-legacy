#
# Pthread
#

if(MSVC)
	set(PTHREAD_LIBRARY_NAME pthreadVCE2)
else(MSVC)
	set(PTHREAD_LIBRARY_NAME pthread)
endif(MSVC)

# Search for Pthreads on user's system.
FIND_LIBRARY(PTHREAD_LIBRARY ${PTHREAD_LIBRARY_NAME} ${MY_PTHREAD_LIB_DIRS})

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PTHREAD DEFAULT_MSG PTHREAD_LIBRARY)
