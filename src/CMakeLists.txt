
find_package(wxWidgets REQUIRED COMPONENTS gl core)
find_package(OpenGL REQUIRED)

include(${wxWidgets_USE_FILE})

add_executable(panoglview
    boundaries.cpp
    boundaries.h
    panoapp.cpp
    panoapp.h
    panocanvas.cpp
    panocanvas.h
    panoframe.cpp
    panoframe.h
    panoglplugin.h
    panointeractivecanvas.cpp
    panointeractivecanvas.h
    panopluginframe.h
)


target_link_libraries(panoglview
    PRIVATE
        ${wxWidgets_LIBRARIES}
        OpenGL::GL
        OpenGL::GLU
)
