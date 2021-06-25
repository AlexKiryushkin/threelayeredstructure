
cmake_minimum_required(VERSION 3.10)
project(tls)

add_library( tls_lib SHARED calculate_unknown_diffusivity.cpp)
target_compile_definitions(tls_lib PUBLIC tls_lib)
target_compile_definitions(tls_lib PRIVATE TLS_EXPORTS)

add_executable( tls main.cpp
                    read_parameters.cpp )
target_compile_features(tls PRIVATE cxx_std_14)
target_link_libraries( tls PUBLIC tls_lib )