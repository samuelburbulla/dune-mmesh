dune_symlink_to_source_files(FILES
  grids
)

dune_add_test(NAME test-mmesh-2d SOURCES test-mmesh-2d.cc)
target_link_libraries("test-mmesh-2d" CGAL)

dune_add_test(NAME test-mmesh-3d SOURCES test-mmesh-3d.cc)
target_link_libraries("test-mmesh-3d" CGAL)

dune_add_test(NAME test-mmesh-structured-2d SOURCES test-mmesh-structured-2d.cc)
target_link_libraries("test-mmesh-structured-2d" CGAL)
