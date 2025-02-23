cmake_minimum_required(VERSION 3.5)

project(hw1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the glm submodule as a dependency
add_subdirectory(glm)

# Add the executable target for your project
add_executable(hw1 main.cpp)

# Link glm to your target
target_link_libraries(hw1 PRIVATE glm::glm)

# Install configuration
include(GNUInstallDirs)
install(TARGETS hw1
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
