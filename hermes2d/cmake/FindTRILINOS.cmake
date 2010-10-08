#
# Trilinos
#
# Looks for Trilinos packages
#
# Required packages are:
# - Epetra, Teuchos
#
SET(MY_TRILINOS_LIB_DIRS $ENV{MY_TRILINOS_LIB_DIRS})
SET(MY_TRILINOS_INC_DIRS $ENV{MY_TRILINOS_INC_DIRS})

IF(NOT MY_TRILINOS_LIB_DIRS OR NOT MY_TRILINOS_INC_DIRS)
  SET(MY_TRILINOS_LIB_DIRS $ENV{SAGE_LOCAL}/lib)
  SET(MY_TRILINOS_INC_DIRS $ENV{SAGE_LOCAL}/include)  
ENDIF(NOT MY_TRILINOS_LIB_DIRS OR NOT MY_TRILINOS_INC_DIRS)

FIND_PATH(AMESOS_INCLUDE_PATH       Amesos.h             ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(AZTECOO_INCLUDE_PATH      AztecOO.h            ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(EPETRA_INCLUDE_PATH       Epetra_Object.h      ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(IFPACK_INCLUDE_PATH       Ifpack.h             ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(LOCA_INCLUDE_PATH         LOCA.H               ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(ML_INCLUDE_PATH           MLAPI.h              ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(NOX_INCLUDE_PATH          NOX.H                ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(TEUCHOS_INCLUDE_PATH      Teuchos_Object.hpp   ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(KOMPLEX_INCLUDE_PATH      Komplex_Version.h    ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)

FIND_PATH(LOCA_EPETRA_INCLUDE_PATH  LOCA_Epetra.H        ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(NOX_EPETRA_INCLUDE_PATH   NOX_Epetra.H         ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)
FIND_PATH(EPETRAEXT_INCLUDE_PATH    EpetraExt_Version.h  ${MY_TRILINOS_INC_DIRS}  NO_DEFAULT_PATH)

FIND_LIBRARY(AMESOS_LIBRARY         amesos               ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(AZTECOO_LIBRARY        aztecoo              ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(EPETRA_LIBRARY         epetra               ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(IFPACK_LIBRARY         ifpack               ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(LOCA_LIBRARY           loca                 ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(ML_LIBRARY             ml                   ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(NOX_LIBRARY            nox                  ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(TEUCHOS_LIBRARY        teuchos              ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(KOMPLEX_LIBRARY        komplex              ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)

FIND_LIBRARY(LOCA_EPETRA_LIBRARY    locaepetra           ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(NOX_EPETRA_LIBRARY     noxepetra            ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)
FIND_LIBRARY(EPETRAEXT_LIBRARY      epetraext            ${MY_TRILINOS_LIB_DIRS}  NO_DEFAULT_PATH)

INCLUDE(FindPackageHandleStandardArgs)

IF(EPETRA_INCLUDE_PATH AND EPETRA_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${EPETRA_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${EPETRA_LIBRARY})
	SET(HAVE_EPETRA YES)
ENDIF(EPETRA_INCLUDE_PATH AND EPETRA_LIBRARY)

IF(TEUCHOS_INCLUDE_PATH AND TEUCHOS_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${TEUCHOS_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${TEUCHOS_LIBRARY})
	SET(HAVE_TEUCHOS YES)
ENDIF(TEUCHOS_INCLUDE_PATH AND TEUCHOS_LIBRARY)

find_package_handle_standard_args(EPETRA DEFAULT_MSG EPETRA_LIBRARY)
find_package_handle_standard_args(TEUCHOS DEFAULT_MSG TEUCHOS_LIBRARY)

IF(AMESOS_INCLUDE_PATH AND AMESOS_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${AMESOS_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${AMESOS_LIBRARY})
	SET(HAVE_AMESOS YES)
	find_package_handle_standard_args(AMESOS DEFAULT_MSG AMESOS_LIBRARY)
ENDIF(AMESOS_INCLUDE_PATH AND AMESOS_LIBRARY)

IF(AZTECOO_INCLUDE_PATH AND AZTECOO_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${AZTECOO_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${AZTECOO_LIBRARY})
	SET(HAVE_AZTECOO YES)
	find_package_handle_standard_args(AZTECOO DEFAULT_MSG AZTECOO_LIBRARY)
ENDIF(AZTECOO_INCLUDE_PATH AND AZTECOO_LIBRARY)

IF(IFPACK_INCLUDE_PATH AND IFPACK_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${IFPACK_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${IFPACK_LIBRARY})
	SET(HAVE_IFPACK YES)
	find_package_handle_standard_args(IFPACK DEFAULT_MSG IFPACK_LIBRARY)
ENDIF(IFPACK_INCLUDE_PATH AND IFPACK_LIBRARY)

IF(LOCA_INCLUDE_PATH AND LOCA_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${LOCA_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${LOCA_LIBRARY})
	SET(HAVE_LOCA YES)
	find_package_handle_standard_args(LOCA DEFAULT_MSG LOCA_LIBRARY)
ENDIF(LOCA_INCLUDE_PATH AND LOCA_LIBRARY)

IF(ML_INCLUDE_PATH AND ML_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${ML_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${ML_LIBRARY})
	SET(HAVE_ML YES)
	find_package_handle_standard_args(ML DEFAULT_MSG ML_LIBRARY)
ENDIF(ML_INCLUDE_PATH AND ML_LIBRARY)

IF(NOX_INCLUDE_PATH AND NOX_LIBRARY AND NOX_EPETRA_INCLUDE_PATH AND NOX_EPETRA_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${NOX_INCLUDE_PATH} ${NOX_EPETRA_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${NOX_LIBRARY} ${NOX_EPETRA_LIBRARY})
	SET(HAVE_NOX YES)
	find_package_handle_standard_args(NOX DEFAULT_MSG NOX_LIBRARY)
ENDIF(NOX_INCLUDE_PATH AND NOX_LIBRARY AND NOX_EPETRA_INCLUDE_PATH AND NOX_EPETRA_LIBRARY)

IF(EPETRAEXT_INCLUDE_PATH AND EPETRAEXT_LIBRARY)
	SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${EPETRAEXT_INCLUDE_PATH})
	SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${EPETRAEXT_LIBRARY})
	SET(HAVE_EPETRAEXT YES)
	find_package_handle_standard_args(EPETRAEXT DEFAULT_MSG EPETRAEXT_LIBRARY)
ENDIF(EPETRAEXT_INCLUDE_PATH AND EPETRAEXT_LIBRARY)

IF(H2D_COMPLEX)
	IF(KOMPLEX_INCLUDE_PATH AND KOMPLEX_LIBRARY)
		SET(TRILINOS_INCLUDE_DIR ${TRILINOS_INCLUDE_DIR} ${KOMPLEX_INCLUDE_PATH})
		SET(TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES} ${KOMPLEX_LIBRARY})
		SET(HAVE_KOMPLEX YES)
		find_package_handle_standard_args(KOMPLEX DEFAULT_MSG KOMPLEX_LIBRARY)
	ENDIF(KOMPLEX_INCLUDE_PATH AND KOMPLEX_LIBRARY)
ENDIF(H2D_COMPLEX)

LIST(REMOVE_DUPLICATES TRILINOS_INCLUDE_DIR)


IF(EPETRA_FOUND AND TEUCHOS_FOUND)
	SET(TRILINOS_FOUND TRUE)
ENDIF(EPETRA_FOUND AND TEUCHOS_FOUND)

IF(H2D_COMPLEX)
	# Komplex is a crucial for compiling complex version of Hermes with Trilinos support.
	IF(NOT KOMPLEX_FOUND)
        MESSAGE(STATUS "Komplex not found.")
		SET(TRILINOS_FOUND FALSE)
	ENDIF(NOT KOMPLEX_FOUND)
ENDIF(H2D_COMPLEX)

IF(TRILINOS_FOUND)
	MESSAGE(STATUS "Trilinos packages found.")
ELSE (TRILINOS_FOUND)
	IF (TRILINOS_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find Trilinos or one of its packages.")
	ENDIF (TRILINOS_FIND_REQUIRED)
ENDIF(TRILINOS_FOUND)
