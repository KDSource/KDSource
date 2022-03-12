cmake_minimum_required(VERSION 3.0.0)

add_executable(test1 "test1.c")
target_link_libraries(test1 PUBLIC kdsource)
add_test(NAME Test1 COMMAND test1)
