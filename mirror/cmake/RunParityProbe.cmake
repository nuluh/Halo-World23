foreach(required_var IN ITEMS PROBE OUTPUT PYTHON COMPARE SOURCE_DIR)
    if(NOT DEFINED ${required_var} OR "${${required_var}}" STREQUAL "")
        message(FATAL_ERROR "RunParityProbe.cmake requires -D${required_var}=...")
    endif()
endforeach()

get_filename_component(output_dir "${OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")

execute_process(
    COMMAND "${PROBE}"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    OUTPUT_FILE "${OUTPUT}"
    RESULT_VARIABLE probe_result
)
if(NOT probe_result EQUAL 0)
    message(FATAL_ERROR "Parity probe failed with exit code ${probe_result}")
endif()

execute_process(
    COMMAND "${PYTHON}" "${COMPARE}" "${OUTPUT}"
    WORKING_DIRECTORY "${SOURCE_DIR}"
    RESULT_VARIABLE compare_result
)
if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "Frozen C++/Fortran comparison failed with exit code ${compare_result}")
endif()
