
add_library(CHSpline CHSpline.cpp)

add_executable (CHSplineDemo main.cpp CHSpline.cpp)
target_link_libraries (CHSplineDemo LINK_PUBLIC CHSpline)
