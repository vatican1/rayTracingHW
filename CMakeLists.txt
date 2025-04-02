cmake_minimum_required(VERSION 3.5)

project(RayTracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the glm submodule as a dependency
add_subdirectory(glm)

find_package(OpenMP REQUIRED)

# Add the executable target for your project
add_executable(RayTracer main.cpp)

# Link glm to your target
target_link_libraries(RayTracer PRIVATE glm::glm)

# Install configuration
include(GNUInstallDirs)
install(TARGETS RayTracer
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
target_link_libraries(RayTracer PRIVATE OpenMP::OpenMP_CXX)
