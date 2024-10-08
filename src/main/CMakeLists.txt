#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2019-2020,2022-2023 Omics Data Automation, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

set(GenomicsDB_library_sources 
    cpp/src/query_operations/variant_operations.cc
    cpp/src/query_operations/broad_combined_gvcf.cc
    cpp/src/query_operations/variant_operations_columnar.cc
    cpp/src/query_operations/broad_combined_gvcf_columnar.cc
    cpp/src/query_operations/gt_remapper_gvcf_columnar_iterator.cc
    cpp/src/genomicsdb/variant_cell.cc
    cpp/src/genomicsdb/variant_storage_manager.cc
    cpp/src/genomicsdb/variant_field_data.cc
    cpp/src/genomicsdb/variant_array_schema.cc
    cpp/src/genomicsdb/variant_field_handler.cc
    cpp/src/genomicsdb/variant.cc
    cpp/src/genomicsdb/query_variants.cc
    cpp/src/genomicsdb/genomicsdb_columnar_field.cc
    cpp/src/genomicsdb/genomicsdb_iterators.cc
    cpp/src/genomicsdb/genomicsdb_multid_vector_field.cc
    cpp/src/loader/tiledb_loader_text_file.cc
    cpp/src/loader/load_operators.cc
    cpp/src/loader/genomicsdb_importer.cc
    cpp/src/loader/tiledb_loader_file_base.cc
    cpp/src/loader/tiledb_loader.cc
    cpp/src/utils/memory_measure.cc
    cpp/src/utils/histogram.cc
    cpp/src/utils/vid_mapper_pb.cc
    cpp/src/utils/lut.cc
    cpp/src/utils/known_field_info.cc
    cpp/src/utils/vid_mapper.cc
    cpp/src/utils/vid_mapper_sql.cc
    cpp/src/utils/timer.cc
    cpp/src/utils/genomicsdb_logger.cc
    cpp/src/utils/alleles_combiner_gvcf_columnar_iterator.cc
    cpp/src/vcf/vcf_adapter.cc
    cpp/src/vcf/genomicsdb_bcf_generator.cc
    cpp/src/vcf/vcf2binary.cc
    cpp/src/vcf/hfile_genomicsdb.c
    cpp/src/vcf/htslib_fs_adapter.cc
    cpp/src/config/variant_query_config.cc
    cpp/src/config/genomicsdb_config_base.cc
    cpp/src/config/json_config.cc
    cpp/src/config/pb_config.cc
    cpp/src/api/annotation_service.cc
    cpp/src/api/genomicsdb.cc
    cpp/src/api/genomicsdb_field.cc
    cpp/src/api/genomicsdb_plink.cc
    cpp/src/api/genomicsdb_utils.cc
    cpp/src/api/genomicsdb_json_processor.cc
    )

if(BUILD_NANOARROW)
  set(GenomicsDB_library_sources ${GenomicsDB_library_sources}
    cpp/src/api/genomicsdb_nanoarrow_processor.cc
    )
endif()    

include_directories(${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS})
set(GenomicsDB_library_sources
    ${GenomicsDB_library_sources}
    ${PROTOBUF_GENERATED_CXX_SRCS}
    )

if(PROTOBUF_REGENERATE)
    #Must be set here - see https://cmake.org/cmake/help/v3.3/command/set_source_files_properties.html
    set_source_files_properties(${PROTOBUF_GENERATED_CXX_SRCS} ${PROTOBUF_GENERATED_CXX_HDRS} ${PROTOBUF_GENERATED_PYTHON_SRCS} ${PROTOBUF_GENERATED_GO_SRCS} ${PROTOBUF_GENERATED_R_SRCS} PROPERTIES GENERATED True)
endif()

if(BUILD_JAVA)
    set(GenomicsDB_library_sources ${GenomicsDB_library_sources}
        jni/src/genomicsdb_GenomicsDBUtils.cc
        jni/src/genomicsdb_GenomicsDBImporter.cc
        jni/src/genomicsdb_GenomicsDBQuery.cc
        jni/src/genomicsdb_GenomicsDBQueryStream.cc
        jni/src/genomicsdb_jni_init.cc
        )
endif()

#Compile sources with PIC
add_library(GenomicsDB_library_object_files OBJECT ${GenomicsDB_library_sources})
if(PROTOBUF_REGENERATE)
    add_dependencies(GenomicsDB_library_object_files PROTOBUF_GENERATED_CXX_TARGET)
endif()
set_property(TARGET GenomicsDB_library_object_files PROPERTY POSITION_INDEPENDENT_CODE ON)

#Create the shared and static libraries
add_library(genomicsdb STATIC $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
#Link in other libraries
add_library(tiledbgenomicsdb SHARED $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
if(TILEDB_SOURCE_DIR)
    target_link_libraries(tiledbgenomicsdb tiledb_static)
else()
    target_link_libraries(tiledbgenomicsdb ${TILEDB_LIBRARY})
endif()
if(HTSLIB_SOURCE_DIR)
  add_dependencies(tiledbgenomicsdb htslib)
endif()
target_link_libraries(tiledbgenomicsdb ${HTSLIB_LIBRARY} ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES})
if(LIBRT_LIBRARY)
    target_link_libraries(tiledbgenomicsdb ${LIBRT_LIBRARY})
endif()
if(BUILD_NANOARROW)
  add_dependencies(tiledbgenomicsdb nanoarrow_ep)
endif()
if(LIBCSV_FOUND)
    target_link_libraries(tiledbgenomicsdb ${LIBCSV_LIBRARY})
endif()
#if(BUILD_DISTRIBUTABLE_LIBRARY)
    #target_link_libraries(tiledbgenomicsdb ${OPENSSL_LIBRARIES})
#endif()

set_target_properties(tiledbgenomicsdb PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})

install(TARGETS genomicsdb tiledbgenomicsdb LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES cpp/include/api/genomicsdb.h  cpp/include/api/genomicsdb_exception.h cpp/include/api/genomicsdb_status.h cpp/include/api/genomicsdb_utils.h
        DESTINATION include)
if(BUILD_FOR_PYTHON)
  install(FILES ${PROTOBUF_GENERATED_PYTHON_SRCS} DESTINATION genomicsdb/protobuf/python)
endif()
if(BUILD_FOR_GO)
  install(FILES ${PROTOBUF_GENERATED_GO_SRCS} DESTINATION genomicsdb/protobuf/go)
endif()
if (BUILD_FOR_R)
  install(FILES ${PROTOBUF_GENERATED_GO_SRCS} DESTINATION genomicsdb/protobuf/r)
endif()

get_filename_component(GENOMICSDB_OPENSSL_ROOT_DIR ${OPENSSL_CRYPTO_LIBRARY}/../ ABSOLUTE)
configure_file(genomicsdb.pc.in genomicsdb.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/genomicsdb.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
