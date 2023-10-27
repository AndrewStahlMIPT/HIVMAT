# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Google_tests) # это строчка не нужна


file(GLOB files "*.cpp" )
foreach(filename ${files})
    get_filename_component (TName ${filename} NAME_WE)
    add_executable("${TName}" ${filename} num_diff_tests.cpp ../SOURCES/picatelonlain.hpp
            integratoror_tests.cpp)
    add_test(NAME "${TName}" COMMAND ${filename})
    target_link_libraries(${TName} HIVMAT gtest gtest_main)
endforeach()