cmake_minimum_required(VERSION 2.8)

project(sfview)

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})
link_directories(${PNG_LIBRARY_DIRS})
add_definitions(${PNG_DEFINITIONS})

add_executable(sfview
  glwindow.cc
  mesh-surface.cc
  nurbs-surface.cc
  sfview.cc
  surface.cc
  utilities.cc
  )

target_link_libraries(sfview ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${PNG_LIBRARIES})
