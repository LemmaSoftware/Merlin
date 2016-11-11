add_executable( KernelV0 KernelV0.cpp  )
target_link_libraries(  KernelV0  "lemmacore" "fdem1d" "merlin")

INSTALL_TARGETS( "/share/Merlin/"
	KernelV0
)
