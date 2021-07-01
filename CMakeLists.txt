
cmake_minimum_required(VERSION 3.10)
project(tls)

MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
  IF(MSVC)
    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(PrecompiledBinary "${CMAKE_CURRENT_BINARY_DIR}/${PrecompiledBasename}.pch")
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_OUTPUTS "${PrecompiledBinary}")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledHeader}\" /FI\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_DEPENDS "${PrecompiledBinary}")  
    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})
  ENDIF(MSVC)
ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)

list(APPEND MySources "calculate_unknown_diffusivity.cpp" "dllmain.cpp")
 

ADD_MSVC_PRECOMPILED_HEADER("pch.h" "pch.cpp" MySources)

add_library( tls_lib SHARED ${MySources})


target_compile_definitions(tls_lib PUBLIC tls_lib)
target_compile_definitions(tls_lib PRIVATE TLS_EXPORTS)

add_executable( tls main.cpp
                    read_parameters.cpp )
target_compile_features(tls PRIVATE cxx_std_14)
target_link_libraries( tls PUBLIC tls_lib )