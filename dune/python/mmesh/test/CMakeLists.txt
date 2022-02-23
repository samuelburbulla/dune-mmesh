dune_python_add_test(NAME interfaceindicator
                     SCRIPT interfaceindicator.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)

dune_python_add_test(NAME normals
                     SCRIPT normals.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)

dune_python_add_test(NAME distance
                     SCRIPT distance.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
