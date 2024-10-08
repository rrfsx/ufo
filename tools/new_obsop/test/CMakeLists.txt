file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/testinput")
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
  "${CMAKE_CURRENT_SOURCE_DIR}/testinput/autogenerated.yaml"
  "${CMAKE_CURRENT_BINARY_DIR}/testinput/autogenerated.yaml")

# Note: can't use filters_testdata.nc4 since it doesn't contain latitudes or longitudes.
# link the local test Data directory into the main ufo test Data directory
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
  "${CMAKE_CURRENT_BINARY_DIR}/../../../test/Data"
  "${CMAKE_CURRENT_BINARY_DIR}/Data")

# Test the generated operator.
ecbuild_add_test(TARGET  test_ufo_opr_autogenerated
                 SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../../../test/mains/TestObsOperator.cc"
                         $<TARGET_OBJECTS:obsautogenerated>
                 ARGS    "testinput/autogenerated.yaml"
                 ENVIRONMENT OOPS_TRAPFPE=1
                 LIBS    ufo )

# (We can't test the linear operator as long as its implementation is empty, since the test
# expects its output to vary as a function of input data.)
