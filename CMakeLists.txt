
cmake_minimum_required(VERSION 3.10)
project(tls)

add_executable( tls main.cpp
                    calculate_unknown_diffusivity.cpp
                    read_parameters.cpp )
target_compile_features(tls PRIVATE cxx_std_14)