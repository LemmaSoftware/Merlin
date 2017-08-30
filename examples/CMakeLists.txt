add_executable( KernelV0 KernelV0.cpp  )
target_link_libraries(  KernelV0  "lemmacore" "fdem1d" "merlin")

add_executable( KernelV0-2 KernelV0-2.cpp  )
target_link_libraries(  KernelV0-2  "lemmacore" "fdem1d" "merlin")

add_executable( KV0-3loops KV0-3loops.cpp  )
target_link_libraries(  KV0-3loops  "lemmacore" "fdem1d" "merlin")

add_executable( ModelAligner ModelAligner.cpp ) 
target_link_libraries(  ModelAligner  "lemmacore" "fdem1d" "merlin")

add_executable( ForwardFID ForwardFID.cpp  )
target_link_libraries(  ForwardFID  "lemmacore" "fdem1d" "merlin")

add_executable( Coupling Coupling.cpp  )
target_link_libraries(  Coupling  "lemmacore" "fdem1d" "merlin")

add_executable( Interference Interference.cpp  )
target_link_libraries(  Interference  "lemmacore" "fdem1d" "merlin")

# Linking
if ( LEMMA_VTK6_SUPPORT OR LEMMA_VTK7_SUPPORT ) 
	target_link_libraries( KernelV0 ${VTK_LIBRARIES})
	target_link_libraries( KV0-3loops ${VTK_LIBRARIES})
	target_link_libraries( Coupling ${VTK_LIBRARIES})
	target_link_libraries( Interference ${VTK_LIBRARIES})
endif()

INSTALL_TARGETS( "/share/Merlin/"
	KernelV0
	Coupling
	Interference
)
