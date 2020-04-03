
if( DUNE_ENABLE_PYTHONBINDINGS )
  dune_python_add_test(NAME test-python-mmesh
                       COMMAND ${PYTHON_EXECUTABLE} mmesh.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(mmesh.py mmesh.py COPYONLY)
endif()
