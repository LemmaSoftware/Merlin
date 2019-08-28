
add_library(pyMerlin MODULE pyMerlin.cpp)
target_link_libraries(pyMerlin PRIVATE pybind11::module lemmacore fdem1d merlin )
set_target_properties(pyMerlin PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                       SUFFIX "${PYTHON_MODULE_EXTENSION}"
                                                      OUTPUT_NAME "Merlin"
)

install(TARGETS pyMerlin
	COMPONENT python
	RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/pyLemma/pyLemma/"
	LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/pyLemma/pyLemma/"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/pyLemma/pyLemma/"
)
