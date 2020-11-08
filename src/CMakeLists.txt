add_executable(sfml_cuda
    main.cpp
    OptixDenoiser.cpp
    options.cpp
    event_handler.cpp
    scenes/full_screen_opengl.cpp

    # scenes/full_screen_opengl.cu
)
#Optix

set(OptiX_INSTALL_DIR "C:/ProgramData/NVIDIA\ Corporation/OptiX\ SDK\ 7.2.0/SDK" CACHE PATH "Path to OptiX installed location.")

# Search for the OptiX libraries and include files.
#find_package(OptiX REQUIRED)

# Add the path to the OptiX headers to our include paths.
include_directories(
  "C:/ProgramData/NVIDIA\ Corporation/OptiX\ SDK\ 7.2.0/include"
  "C:/ProgramData/NVIDIA\ Corporation/OptiX\ SDK\ 7.2.0/SDK"
  "${OptiX_INSTALL_DIR}/cuda"
  )


#enable profiler
set(ENABLE_OPTICK OFF)

include_directories( ${CMAKE_SOURCE_DIR}/include )
if (ENABLE_OPTICK)
	add_definitions(-DUSE_OPTICK)
  if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
    target_link_libraries(sfml_cuda PRIVATE ${CMAKE_SOURCE_DIR}/lib/x64/debug/OptickCore.lib)
    add_custom_command(
        TARGET sfml_cuda POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_SOURCE_DIR}/lib/x64/debug/OptickCore.dll
          $<TARGET_FILE_DIR:sfml_cuda>)

  else()
    target_link_libraries(sfml_cuda PRIVATE ${CMAKE_SOURCE_DIR}/lib/x64/release/OptickCore.lib)
    add_custom_command(
        TARGET sfml_cuda POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_SOURCE_DIR}/lib/x64/release/OptickCore.dll
          $<TARGET_FILE_DIR:sfml_cuda>)

  endif()

	set(COMMON_LIBS ${COMMON_LIBS} ProfilerCore64)
endif()

add_custom_command(TARGET sfml_cuda POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_SOURCE_DIR}/media $<TARGET_FILE_DIR:sfml_cuda>)


target_compile_options(sfml_cuda PRIVATE                        
    $<$<COMPILE_LANGUAGE:CUDA>:${ARCH_FLAGS} --expt-relaxed-constexpr>
)
set_target_properties(sfml_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_include_directories(sfml_cuda SYSTEM PUBLIC ${CUDA_INCLUDE_DIRS})
target_link_libraries(sfml_cuda PRIVATE project_options project_warnings ${CUDA_LIBRARIES})

target_compile_options(sfml_cuda PUBLIC /Zc:twoPhase-) 

target_link_libraries_system(
  sfml_cuda PRIVATE  
  CONAN_PKG::docopt.cpp
  CONAN_PKG::fmt
  CONAN_PKG::spdlog
  CONAN_PKG::imgui-sfml
  CONAN_PKG::glew
  CONAN_PKG::glad
  CONAN_PKG::eigen
)

