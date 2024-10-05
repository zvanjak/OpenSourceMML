cmake_minimum_required(VERSION 3.20)

set (CMAKE_CXX_STANDARD 20)

project(MML_DemoApp VERSION 1.0)

set(MML_DEMO_APP_SOURCES  demo_app_main.cpp 
                          demo_app/demo_vector.cpp
)

add_executable(MML_DemoApp ${MML_DEMO_APP_SOURCES})
target_include_directories(MML_DemoApp PRIVATE ../include)


project(MML_DocsApp VERSION 1.0)
set(MML_DOCS_APP_SOURCES    docs_demos_main.cpp

                            docs_demos/docs_demo_vector.cpp
)
add_executable(MML_DocsApp ${MML_DOCS_APP_SOURCES} )
target_include_directories(MML_DocsApp PRIVATE ../include)


project(MML_ReadmeApp VERSION 1.0)
set(MML_README_APP_SOURCES  readme_examples_main.cpp

                            readme_examples/readme1_vectors_matrices.cpp
)
add_executable(MML_ReadmeApp ${MML_README_APP_SOURCES} )
target_include_directories(MML_ReadmeApp PRIVATE ../include)


project(MML_ExamplesApp VERSION 1.0)
set(MML_EXAMPLES_APP_SOURCES  examples_main.cpp

                              examples/example1_kosi_hitac.cpp
)
add_executable(MML_ExamplesApp ${MML_EXAMPLES_APP_SOURCES} )
target_include_directories(MML_ExamplesApp PRIVATE ../include)
