# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/include/htslib"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src/htslib-build"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/tmp"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src/htslib-stamp"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src"
  "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src/htslib-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src/htslib-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/chaokuan-hao/Documents/Projects/Junction-Extractor/build/htslib-prefix/src/htslib-stamp${cfgdir}") # cfgdir has leading slash
endif()
