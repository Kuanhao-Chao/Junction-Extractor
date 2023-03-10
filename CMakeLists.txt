cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(junction_extractor)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -DNOCURL=1")

include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)

ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --without-libdeflate --disable-libcurl --disable-lzma
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(junction_extractor src/junction_extractor.cpp 
        src/GSam.cpp
        src/GSam.h
        src/tmerge.cpp
        src/tmerge.h
        # src/commons.h
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gff.h

        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/gdna.h

        ${CMAKE_SOURCE_DIR}/include/gclib/codons.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/codons.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GFaSeqGet.h

        ${CMAKE_SOURCE_DIR}/include/gclib/GList.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GHashMap.hh
        ${CMAKE_SOURCE_DIR}/include/gclib/GBitVec.h
        
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp
        )

add_dependencies(junction_extractor htslib)
target_link_libraries(junction_extractor ${CMAKE_SOURCE_DIR}/include/htslib/libhts.a z lzma bz2 pthread)
        
set_property(TARGET junction_extractor PROPERTY CXX_STANDARD 17)

install(TARGETS junction_extractor DESTINATION ${CMAKE_INSTALL_PREFIX})