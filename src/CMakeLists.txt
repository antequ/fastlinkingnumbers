set(MODEL_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/curve.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/model.cpp)

set(VERIFY_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/compute_link.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/discretize.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/linking_number_certificate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/potential_link_search.cpp)

project(verifycurves
    VERSION 1.0.0
    DESCRIPTION
    "Topology Verification Tool Using Fast Linking Numbers")

add_executable(verifycurves
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${MODEL_FILES}
    ${VERIFY_FILES})

target_link_libraries(verifycurves
    Eigen3::Eigen
    gflags
    OpenMP::OpenMP_CXX)

project(obj2bcc
    VERSION 1.0.0
    DESCRIPTION
    "Converts Wavefront OBJ files to Binary Curve Collection (BCC) files.")

add_executable(obj2bcc
    ${CMAKE_CURRENT_SOURCE_DIR}/obj2bcc.cpp
    ${MODEL_FILES})

target_link_libraries(obj2bcc
    Eigen3::Eigen
    gflags)

project(replicability
    VERSION 1.0.0
    DESCRIPTION
    "Code Replicability Demo for Generating Table 1 in Paper")

add_executable(replicability
    ${CMAKE_CURRENT_SOURCE_DIR}/replicability_script.cpp
    ${MODEL_FILES}
    ${VERIFY_FILES})

target_link_libraries(replicability
    Eigen3::Eigen
    gflags
    OpenMP::OpenMP_CXX)

project(libverifycurves
    VERSION 1.0.0
    DESCRIPTION
    "Topology Verification Library Using Fast Linking Numbers")

add_library(libverifycurves SHARED
    ${MODEL_FILES}
    ${VERIFY_FILES})

target_link_libraries(libverifycurves
    Eigen3::Eigen
    OpenMP::OpenMP_CXX)

set_target_properties(libverifycurves
    PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(libverifycurves
    PROPERTIES OUTPUT_NAME verifycurves)

