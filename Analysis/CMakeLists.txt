set(SOURCES
    src/functions.cc
)

set(HEADERS
    include/functions.h
)

add_library(Analysis SHARED ${SOURCES} ${HEADERS})

target_include_directories(Analysis PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(Analysis PUBLIC 
    MCPSim
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Physics
    ROOT::MathCore
)

set_target_properties(Analysis
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

add_executable(Analysis_exe Analysis.cc)
target_link_libraries(Analysis_exe PRIVATE Analysis MCPSim)
set_target_properties(Analysis_exe
    PROPERTIES
    OUTPUT_NAME "Analysis"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties(Analysis_exe PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)

install(TARGETS Analysis Analysis_exe
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include
)

install(FILES ${HEADERS}
    DESTINATION include/Analysis
)
