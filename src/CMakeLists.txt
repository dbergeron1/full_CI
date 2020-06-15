cmake_minimum_required (VERSION 2.6)

project(full_CI)

add_executable(full_CI graph_2D.cpp full_CI_main.cpp)

target_compile_options(full_CI PRIVATE -std=c++11 -O0)

target_include_directories(full_CI PRIVATE /usr/local/opt/llvm/include /usr/local/Cellar/llvm/7.0.0/lib/clang/7.0.0/include)

target_link_libraries(full_CI gmpxx gmp blas lapack)

install(TARGETS full_CI DESTINATION bin)

