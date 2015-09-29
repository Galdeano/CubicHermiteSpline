
add_library(CHSpline SHARED CHSpline.cpp)

add_executable (CHSplineDemo main.cpp CHSpline.cpp)
target_link_libraries (CHSplineDemo LINK_PUBLIC CHSpline)

install(TARGETS CHSplineDemo 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
