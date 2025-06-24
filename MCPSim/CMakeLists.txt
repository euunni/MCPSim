file(GLOB SOURCES src/*.cc)
file(GLOB HEADERS include/*.h)

# Set Eigen path
set(EIGEN3_INCLUDE_DIR "/home/jangh/eigen-3.4.0")

find_package(ROOT REQUIRED)

# ROOT 딕셔너리 생성
ROOT_GENERATE_DICTIONARY(MCPEventDict
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/MCPEventData.h
                         LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/MCPEventDataLinkDef.h
                         MODULE MCPSim
                         OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/include
                         NOINSTALL)

list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/MCPEventDict.cxx)

add_library(MCPSim SHARED ${SOURCES})

target_include_directories(MCPSim PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${EIGEN3_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(MCPSim PUBLIC ${ROOT_LIBRARIES})

set(MCPSIM_CONFIG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/config)

target_compile_definitions(MCPSim PUBLIC
    MCPSIM_CONFIG_DIR="${MCPSIM_CONFIG_DIR}"
)

set_target_properties(MCPSim
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

add_executable(MCPSim_exe MCPSim.cc)
target_link_libraries(MCPSim_exe PRIVATE MCPSim)
set_target_properties(MCPSim_exe
    PROPERTIES
    OUTPUT_NAME "MCPSim"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

set_target_properties(MCPSim_exe PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config 
     DESTINATION ${CMAKE_BINARY_DIR}/MCPSim)

# 설치
install(TARGETS MCPSim MCPSim_exe
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include
)

install(FILES ${HEADERS}
    DESTINATION include/MCPSim
)

# 딕셔너리 파일 설치
install(FILES
    ${CMAKE_BINARY_DIR}/lib/libMCPSim_rdict.pcm
    ${CMAKE_BINARY_DIR}/lib/libMCPSim.rootmap
    DESTINATION lib
)

install(DIRECTORY config
    DESTINATION share/MCPSim)
