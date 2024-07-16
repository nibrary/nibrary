include("${CMAKE_CURRENT_LIST_DIR}/utils.cmake")

conditional_make_directory("${CMAKE_INSTALL_PREFIX}/include/${nibrary}/geogram")
conditional_make_directory("${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram")
conditional_copy_directory("${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram" "${CMAKE_INSTALL_PREFIX}/include/${nibrary}/geogram")
conditional_remove_directory("${CMAKE_INSTALL_PREFIX}/include/geogram1")
