# Configure Merlin 
set(MERLIN_VERSION_MAJOR "0")
set(MERLIN_VERSION_MINOR "1")
set(MERLIN_VERSION_PATCH "0")
set(MERLIN_VERSION "\"${MERLIN_VERSION_MAJOR}.${MERLIN_VERSION_MINOR}.${MERLIN_VERSION_PATCH}\"")
set(MERLIN_VERSION_NOQUOTES "${MERLIN_VERSION_MAJOR}.${MERLIN_VERSION_MINOR}.${MERLIN_VERSION_PATCH}")

configure_file (
	"${CMAKE_CURRENT_SOURCE_DIR}/config/MerlinConfig.h.in"
	"${PROJECT_BINARY_DIR}/include/MerlinConfig.h"
)
install ( FILES  ${PROJECT_BINARY_DIR}/include/MerlinConfig.h   DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Lemma/ )

add_subdirectory("src")
add_library( merlin ${MERLINSOURCE} )  
target_include_directories( merlin PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" )

set_target_properties(merlin PROPERTIES 
	VERSION  "${MERLIN_VERSION_NOQUOTES}"
	SOVERSION "${MERLIN_VERSION_MAJOR}.${MERLIN_VERSION_MINOR}"
	PROJECT_LABEL "MERLIN ${LABEL_SUFFIX}"
	CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

# Linking
target_link_libraries(merlin "lemmacore" "fdem1d" )

# Linking
if ( LEMMA_VTK6_SUPPORT OR LEMMA_VTK7_SUPPORT ) 
	target_link_libraries(merlin ${VTK_LIBRARIES})
endif()

# Testing
if (LEMMA_ENABLE_TESTING)
	add_subdirectory(testing)
endif()

# Python Bindings
if (LEMMA_PYTHON3_BINDINGS)
	add_subdirectory(python)
endif()

# Install
install ( TARGETS merlin DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
install ( FILES include/Merlin  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Lemma ) 
install ( DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Lemma  FILES_MATCHING PATTERN "*.h")

#install ( DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Lemma/  FILES_MATCHING PATTERN "FDEM1D")
#install ( DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/Lemma/FDEM1D  FILES_MATCHING PATTERN "*.h")

# Examples
if (LEMMA_BUILD_EXAMPLES)
	add_subdirectory(examples)
endif()
