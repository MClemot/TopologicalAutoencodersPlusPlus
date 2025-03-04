# overview
function(ttk_print_summary)
    message(STATUS "ttk options -----------------------------------------------------------------")
    message(STATUS "TTK_ENABLE_64BIT_IDS: ${TTK_ENABLE_64BIT_IDS}")
    message(STATUS "TTK_ENABLE_CGAL: ${TTK_ENABLE_CGAL}")
    message(STATUS "TTK_ENABLE_CPU_OPTIMIZATION: ${TTK_ENABLE_CPU_OPTIMIZATION}")
    message(STATUS "TTK_ENABLE_DOUBLE_TEMPLATING: ${TTK_ENABLE_DOUBLE_TEMPLATING}")
    message(STATUS "TTK_ENABLE_EIGEN: ${TTK_ENABLE_EIGEN}")
    message(STATUS "TTK_ENABLE_EMBREE: ${TTK_ENABLE_EMBREE}")
    message(STATUS "TTK_ENABLE_GRAPHVIZ: ${TTK_ENABLE_GRAPHVIZ}")
    message(STATUS "TTK_ENABLE_KAMIKAZE: ${TTK_ENABLE_KAMIKAZE}")
    message(STATUS "TTK_ENABLE_MPI: ${TTK_ENABLE_MPI}")
    message(STATUS "TTK_ENABLE_OPENMP: ${TTK_ENABLE_OPENMP}")
    message(STATUS "TTK_ENABLE_QHULL: ${TTK_ENABLE_QHULL}")
    message(STATUS "TTK_ENABLE_SCIKIT_LEARN: ${TTK_ENABLE_SCIKIT_LEARN}")
    message(STATUS "TTK_ENABLE_SHARED_BASE_LIBRARIES: ${TTK_ENABLE_SHARED_BASE_LIBRARIES}")
    message(STATUS "TTK_ENABLE_SPECTRA: ${TTK_ENABLE_SPECTRA}")
    message(STATUS "TTK_ENABLE_SQLITE3: ${TTK_ENABLE_SQLITE3}")
    message(STATUS "TTK_ENABLE_TORCH: ${TTK_ENABLE_TORCH}")
    message(STATUS "TTK_ENABLE_WEBSOCKETPP: ${TTK_ENABLE_WEBSOCKETPP}")
    message(STATUS "TTK_ENABLE_ZFP: ${TTK_ENABLE_ZFP}")
    message(STATUS "TTK_ENABLE_ZLIB: ${TTK_ENABLE_ZLIB}")
    message(STATUS "ttk build -------------------------------------------------------------------")
    message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
    message(STATUS "TTK_BUILD_DOCUMENTATION: ${TTK_BUILD_DOCUMENTATION}")
    if(TTK_BUILD_DOCUMENTATION)
        message(STATUS "  DOXYGEN_EXECUTABLE: ${DOXYGEN_EXECUTABLE}")
    endif()
    message(STATUS "TTK_BUILD_PARAVIEW_PLUGINS: ${TTK_BUILD_PARAVIEW_PLUGINS}")
    if(TTK_BUILD_PARAVIEW_PLUGINS)
        message(STATUS "  ParaView_DIR: ${ParaView_DIR}")
    endif()
    message(STATUS "TTK_BUILD_STANDALONE_APPS: ${TTK_BUILD_STANDALONE_APPS}")
    message(STATUS "TTK_BUILD_VTK_WRAPPERS: ${TTK_BUILD_VTK_WRAPPERS}")
    if(TTK_BUILD_VTK_WRAPPERS)
        message(STATUS "  VTK_DIR: ${VTK_DIR}")
        message(STATUS "  TTK_BUILD_VTK_PYTHON_MODULE: ${TTK_BUILD_VTK_PYTHON_MODULE}")
    endif()
    message(STATUS "ttk install -----------------------------------------------------------------")
    message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "TTK_PYTHON_MODULE_DIR: ${TTK_PYTHON_MODULE_DIR}")
    message(STATUS "TTK_INSTALL_PLUGIN_DIR: ${TTK_INSTALL_PLUGIN_DIR}")
    message(STATUS "-----------------------------------------------------------------------------")
endfunction()

