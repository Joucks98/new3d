cmake_minimum_required(VERSION 2.8)
PROJECT(New3D)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) 
#find_package(QT5 REQUIRED)
#include(${QT_USE_FILE})

set(Proj_SRCS 
	main.cpp
	new3d.h
	new3d.cpp
	ImageCallBack.h
	ImageCallBack.cpp
	InteractorStylePointCloud.h
	InteractorStylePointCloud.cpp
)
set(Proj_UIS
	new3d.ui
)
set(Proj_Moc_HDRS
	new3d.h
)
qt5_wrap_ui(Proj_UIS_H ${Proj_UIS}
)
qt5_wrap_cpp(Proj_Moc_SRCS ${Proj_Moc_HDRS}
)

include_directories(${Proj_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)

set(Proj_Rcc_SRCS )
add_executable(New3D MACOSX_BUNDLE ${Proj_SRCS} ${Proj_UIS_H} ${Proj_Moc_SRCS} ${Proj_Rcc_SRCS} ) 
target_link_libraries(New3D ${VTK_LIBRARIES} ${QT_LIBRARIES})