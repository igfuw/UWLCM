if(UWLCM_DISABLE)
  message(WARNING "Compilation of some of the options was disabled via the UWLCM_DISABLE flag. Unit tests are most likely to fail.") 
endif()

add_executable(api_test api_test.cpp)
target_compile_features(api_test PRIVATE cxx_std_11)

add_test(api_test_iles api_test ${CMAKE_BINARY_DIR} 1 iles)
add_test(api_test_smg  api_test ${CMAKE_BINARY_DIR} 0 smg " --sgs='smg' ")
add_test(api_test_smgani  api_test ${CMAKE_BINARY_DIR} 0 smgani " --sgs='smgani' ")

# reference data decompression
add_test(NAME SetupReferenceData
         COMMAND tar --zstd -xf ${CMAKE_CURRENT_SOURCE_DIR}/reference_data.tar.zst
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(SetupReferenceData PROPERTIES FIXTURES_SETUP "RefDataFixture")

# Helper to add an "api_test_<scheme>_diff" style test.
function(add_api_test_diff name)
  add_test(api_test_${name}_diff bash -c "
    cd ${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/ && find ./* -maxdepth 0 -type d -print0 |
      while IFS= read -rd '' dir; do
        echo \"\\n============== working with $dir ===============\";
        params=\"$(grep $(echo \"$dir\" | sed 's/..//') ../hash_dict_${name}.txt)\"
        echo \"$params\"
        echo   'comparing const.h5'
        h5diff -v2 --relative=2e-4 --exclude-path \"/misc\" --exclude-path \"/git_revisions\" --exclude-path \"/MPI details\" \"${CMAKE_CURRENT_BINARY_DIR}/output_${name}/$dir/const.h5\" \"${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/$dir/const.h5\"                 &&
        echo   'comparing temp.xmf'                                                                                                                                                      &&
        diff \"${CMAKE_CURRENT_BINARY_DIR}/output_${name}/$dir/temp.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/$dir/temp.xmf\"                                                         &&
        echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
        h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output_${name}/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/$dir/timestep0000000000.h5\" rv     &&
        h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output_${name}/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/$dir/timestep0000000000.h5\" th     &&
        echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
        diff \"${CMAKE_CURRENT_BINARY_DIR}/output_${name}/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_${name}/$dir/timestep0000000000.xmf\"                             || exit 1;
      done
  ")
endfunction()

add_api_test_diff(iles)
add_api_test_diff(smg)
add_api_test_diff(smgani)

# Mark that these tests need the decompressed data
set_tests_properties(api_test_iles_diff api_test_smg_diff api_test_smgani_diff
                     PROPERTIES FIXTURES_REQUIRED "RefDataFixture")


# 3. Define the Cleanup Test (Delete the extracted folder)
add_test(NAME CleanupReferenceData 
         COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_BINARY_DIR}/refdata_iles ${CMAKE_CURRENT_BINARY_DIR}/refdata_smgs ${CMAKE_CURRENT_BINARY_DIR}/refdata_smgani)

# Mark this test as the cleanup step for the fixture
set_tests_properties(CleanupReferenceData PROPERTIES FIXTURES_CLEANUP "RefDataFixture")
