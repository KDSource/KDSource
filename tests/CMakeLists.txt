cmake_minimum_required(VERSION 2.6)

add_executable(test1 "test1.c")
target_link_libraries(test1 PUBLIC ksource)
add_test(NAME Test1 COMMAND test1)
