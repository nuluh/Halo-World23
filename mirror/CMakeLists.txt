cmake_minimum_required(VERSION 3.20)

project(openVFIFE LANGUAGES CXX)

include(CTest)
include(CheckCXXCompilerFlag)

option(OPENVFIFE_ENABLE_HARDENING "Enable project-local native hardening flags" ON)
option(OPENVFIFE_ENABLE_SANITIZERS "Build project validation with ASan and UBSan" OFF)
option(OPENVFIFE_WARNINGS_AS_ERRORS "Treat warnings as errors for openVFIFE-owned targets" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(OPENVFIFE_ENABLE_SANITIZERS AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR
        "OPENVFIFE_ENABLE_SANITIZERS requires a GNU or Clang-family compiler "
        "with AddressSanitizer and UndefinedBehaviorSanitizer support."
    )
endif()

# R1 keeps the bundled Eigen 3.4.0 tree as the offline dependency baseline.
# #47 will own final install/export mechanics; this interface target establishes
# dependency ownership now without downloading or discovering packages at configure time.
add_library(openvfife_eigen INTERFACE)
add_library(openVFIFE::eigen ALIAS openvfife_eigen)
target_include_directories(openvfife_eigen INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/headers>
)

set(OPENVFIFE_CORE_SOURCES
    src/core/particle.cpp
    src/numerics/quadrature.cpp
    src/material/material.cpp
    src/section/section.cpp
    src/section/grid_section.cpp
    src/element/element.cpp
    src/solver/struct_system.cpp
)

add_library(openvfife STATIC ${OPENVFIFE_CORE_SOURCES})
add_library(openVFIFE::openvfife ALIAS openvfife)

target_compile_features(openvfife PUBLIC cxx_std_20)
target_include_directories(openvfife PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(openvfife PUBLIC openVFIFE::eigen)

function(openvfife_add_supported_compile_flag target flag variable_name)
    check_cxx_compiler_flag("${flag}" "${variable_name}")
    if(${variable_name})
        target_compile_options(${target} PRIVATE "${flag}")
    endif()
endfunction()

function(openvfife_apply_project_compile_options target)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(OPENVFIFE_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()

        if(OPENVFIFE_ENABLE_SANITIZERS)
            target_compile_options(${target} PRIVATE
                -O1
                -g
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
            )
        elseif(OPENVFIFE_ENABLE_HARDENING AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
            target_compile_definitions(${target} PRIVATE
                $<$<NOT:$<CONFIG:Debug>>:_FORTIFY_SOURCE=3>
            )
            openvfife_add_supported_compile_flag(${target} -fstack-protector-strong OPENVFIFE_HAVE_STACK_PROTECTOR_STRONG)
            openvfife_add_supported_compile_flag(${target} -fstack-clash-protection OPENVFIFE_HAVE_STACK_CLASH_PROTECTION)
            openvfife_add_supported_compile_flag(${target} -ftrivial-auto-var-init=zero OPENVFIFE_HAVE_TRIVIAL_AUTO_VAR_INIT_ZERO)
        endif()
    elseif(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(OPENVFIFE_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    endif()
endfunction()

function(openvfife_apply_sanitizer_link_options target)
    if(OPENVFIFE_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Keep sanitizer linkage private to openVFIFE validation executables.
        # It is intentionally not part of the public library usage requirements.
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

openvfife_apply_project_compile_options(openvfife)

if(BUILD_TESTING)
    find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)

    function(openvfife_add_test_executable target source)
        add_executable(${target} ${source})
        target_link_libraries(${target} PRIVATE openVFIFE::openvfife)
        target_compile_features(${target} PRIVATE cxx_std_20)
        openvfife_apply_project_compile_options(${target})
        openvfife_apply_sanitizer_link_options(${target})
    endfunction()

    # Each intended public header must compile as the only openVFIFE include in its TU.
    set(OPENVFIFE_PUBLIC_HEADERS
        openvfife/core/types.hpp
        openvfife/core/particle.hpp
        openvfife/numerics/quadrature.hpp
        openvfife/material/uniaxial.hpp
        openvfife/material/experimental.hpp
        openvfife/section/section.hpp
        openvfife/section/experimental.hpp
        openvfife/section/grid_section.hpp
        openvfife/element/base.hpp
        openvfife/element/frame.hpp
        openvfife/element/experimental.hpp
        openvfife/solver/struct_system.hpp
        openvfife/openvfife.hpp
    )
    set(OPENVFIFE_HEADER_SMOKE_TARGETS)
    foreach(public_header IN LISTS OPENVFIFE_PUBLIC_HEADERS)
        string(REPLACE "/" "_" header_target "${public_header}")
        string(REPLACE "." "_" header_target "${header_target}")
        set(header_source "${CMAKE_CURRENT_BINARY_DIR}/header_smoke_${header_target}.cpp")
        file(WRITE "${header_source}" "#include <${public_header}>\n")
        add_library(header_smoke_${header_target} OBJECT "${header_source}")
        target_link_libraries(header_smoke_${header_target} PRIVATE openVFIFE::openvfife)
        target_compile_features(header_smoke_${header_target} PRIVATE cxx_std_20)
        openvfife_apply_project_compile_options(header_smoke_${header_target})
        list(APPEND OPENVFIFE_HEADER_SMOKE_TARGETS header_smoke_${header_target})
    endforeach()
    add_custom_target(openvfife_public_header_smoke DEPENDS ${OPENVFIFE_HEADER_SMOKE_TARGETS})

    function(openvfife_register_native_test target)
        add_test(NAME ${target} COMMAND ${target})
        set_tests_properties(${target} PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            LABELS "native"
        )
    endfunction()

    openvfife_add_test_executable(test_beam2d tests/test_beam2d.cpp)
    openvfife_add_test_executable(test_beam3d tests/test_beam3d.cpp)
    openvfife_add_test_executable(test_element tests/test_element.cpp)
    openvfife_add_test_executable(test_gridsection tests/test_gridsection.cpp)
    openvfife_add_test_executable(test_integrator tests/test_integrator.cpp)
    openvfife_add_test_executable(test_material tests/test_material.cpp)
    openvfife_add_test_executable(test_particle tests/test_particle.cpp)
    openvfife_add_test_executable(test_section tests/test_section.cpp)
    openvfife_add_test_executable(test_structsystem tests/test_structsystem.cpp)

    # The native shim remains at its legacy path until #48, but it now consumes
    # the canonical core target instead of recompiling the seven core sources.
    add_library(openvfife_native_shim OBJECT wasm/vfife_shim.cpp)
    target_link_libraries(openvfife_native_shim PRIVATE openVFIFE::openvfife)
    target_compile_features(openvfife_native_shim PRIVATE cxx_std_20)
    openvfife_apply_project_compile_options(openvfife_native_shim)

    add_executable(test_shim
        tests/test_shim.cpp
        $<TARGET_OBJECTS:openvfife_native_shim>
    )
    target_link_libraries(test_shim PRIVATE openVFIFE::openvfife)
    target_compile_features(test_shim PRIVATE cxx_std_20)
    openvfife_apply_project_compile_options(test_shim)
    openvfife_apply_sanitizer_link_options(test_shim)

    add_executable(test_public_consumer
        tests/consumer/public_consumer_main.cpp
        tests/consumer/public_consumer_model.cpp
    )
    target_link_libraries(test_public_consumer PRIVATE openVFIFE::openvfife)
    target_compile_features(test_public_consumer PRIVATE cxx_std_20)
    openvfife_apply_project_compile_options(test_public_consumer)
    openvfife_apply_sanitizer_link_options(test_public_consumer)
    openvfife_register_native_test(test_public_consumer)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        foreach(surface_header IN ITEMS
            openvfife/openvfife.hpp
            openvfife/solver/struct_system.hpp
        )
            string(REPLACE "/" "_" surface_name "${surface_header}")
            string(REPLACE "." "_" surface_name "${surface_name}")
            add_test(NAME include_surface_${surface_name}
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_public_include_surface.py
                    --compiler ${CMAKE_CXX_COMPILER}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --header ${surface_header}
            )
            set_tests_properties(include_surface_${surface_name} PROPERTIES LABELS "native")
        endforeach()
    endif()

    openvfife_add_test_executable(fortran_parity_probe tests/fortran_parity_probe.cpp)
    openvfife_add_test_executable(frame_progression_probe tests/frame_progression_probe.cpp)

    foreach(test_target IN ITEMS
        test_beam2d
        test_beam3d
        test_element
        test_gridsection
        test_integrator
        test_material
        test_particle
        test_section
        test_structsystem
        test_shim
    )
        openvfife_register_native_test(${test_target})
    endforeach()

    set(OPENVFIFE_PARITY_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/validation")

    add_test(NAME fortran_parity
        COMMAND ${CMAKE_COMMAND}
            -DPROBE=$<TARGET_FILE:fortran_parity_probe>
            -DOUTPUT=${OPENVFIFE_PARITY_OUTPUT_DIR}/fortran_parity_probe.csv
            -DPYTHON=${Python3_EXECUTABLE}
            -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/validation/fortran/golden/compare_cpp.py
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunParityProbe.cmake
    )
    set_tests_properties(fortran_parity PROPERTIES LABELS "native;sanitizer")

    add_test(NAME frame_progression_parity
        COMMAND ${CMAKE_COMMAND}
            -DPROBE=$<TARGET_FILE:frame_progression_probe>
            -DOUTPUT=${OPENVFIFE_PARITY_OUTPUT_DIR}/frame_progression_probe.csv
            -DPYTHON=${Python3_EXECUTABLE}
            -DCOMPARE=${CMAKE_CURRENT_SOURCE_DIR}/validation/fortran/golden/compare_cpp.py
            -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunParityProbe.cmake
    )
    set_tests_properties(frame_progression_parity PROPERTIES LABELS "native;sanitizer")

    # Preserve the established fast sanitizer surface. Beam3D and the shim remain
    # part of normal native validation but were not part of the pre-R1 sanitizer script.
    foreach(sanitizer_test IN ITEMS
        test_beam2d
        test_element
        test_gridsection
        test_integrator
        test_material
        test_particle
        test_section
        test_structsystem
    )
        set_property(TEST ${sanitizer_test} APPEND PROPERTY LABELS sanitizer)
    endforeach()

    add_custom_target(openvfife_sanitizer_targets DEPENDS
        openvfife
        test_beam2d
        test_element
        test_gridsection
        test_integrator
        test_material
        test_particle
        test_section
        test_structsystem
        fortran_parity_probe
        frame_progression_probe
    )
endif()
