add_executable( KernelV0 KernelV0.cpp  )
target_link_libraries(  KernelV0  "lemmacore" "fdem1d" "merlin")

# Linking
if ( LEMMA_VTK6_SUPPORT OR LEMMA_VTK7_SUPPORT ) 
	target_link_libraries( KernelV0 ${VTK_LIBRARIES})
#	target_link_libraries(lemmacore "matplot")
endif()

INSTALL_TARGETS( "/share/Merlin/"
	KernelV0
)
