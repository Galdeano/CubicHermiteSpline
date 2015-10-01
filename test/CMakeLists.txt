
find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

include_directories (../src ${TEST_SOURCE_DIR}/src ${Boost_INCLUDE_DIR}) 

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

enable_testing()

add_executable (test_CHSpline test_CHSpline.cpp ../src/CHSpline.cpp) 

target_link_libraries(test_CHSpline
                      ${Boost_FILESYSTEM_LIBRARY}
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                      gcov
                      )

add_test(MyTest_CHSpline test_CHSpline)
