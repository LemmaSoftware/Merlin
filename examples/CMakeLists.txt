add_executable( KernelV0 KernelV0.cpp  )
target_link_libraries(  KernelV0  "lemmacore" "fdem1d" "merlin")

add_executable( Coupling Coupling.cpp  )
target_link_libraries(  Coupling  "lemmacore" "fdem1d" "merlin")

# Linking
if ( LEMMA_VTK6_SUPPORT OR LEMMA_VTK7_SUPPORT ) 
	target_link_libraries( KernelV0 ${VTK_LIBRARIES})
	target_link_libraries( Coupling ${VTK_LIBRARIES})
endif()

INSTALL_TARGETS( "/share/Merlin/"
	KernelV0
	Coupling
)
