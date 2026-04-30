if(UWLCM_DISABLE)
  message(WARNING "Compilation of some of the options was disabled via the UWLCM_DISABLE flag. Unit tests are most likely to fail.") 
endif()

add_executable(api_test api_test.cpp)
target_compile_features(api_test PRIVATE cxx_std_11)

add_test(api_test_iles api_test ${CMAKE_BINARY_DIR} 1)
add_test(api_test_smg  api_test ${CMAKE_BINARY_DIR} 0 " --sgs=1 ")

# reference data decompression
add_test(NAME SetupReferenceData
         COMMAND tar --zstd -xf ${CMAKE_CURRENT_SOURCE_DIR}/reference_data.tar.zst
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(SetupReferenceData PROPERTIES FIXTURES_SETUP "RefDataFixture")


add_test(api_test_iles_diff bash -c " 
  cd ${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/ && find ./* -maxdepth 0 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"\\n============== working with $dir ===============\";
      params=$(grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt)
      echo \"$params\"
      
      # useful for skipping some tests
      #if [[ $params == *\"cumulus\"* ]]; then
      #  echo \"skipping\"
      #  continue
      #fi

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=2e-4 --exclude-path \"/misc\" --exclude-path \"/git_revisions\" --exclude-path \"/MPI details\"  \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"  \"${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/$dir/const.h5\"   &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/$dir/temp.xmf\"                                                         && 
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/$dir/timestep0000000000.h5\" rv     &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/$dir/timestep0000000000.h5\" th     &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_iles/$dir/timestep0000000000.xmf\"                             || exit 1; 
    done  
")

add_test(api_test_smg_diff bash -c " 
  cd ${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/ && find ./* -maxdepth 0 -type d -print0 | 
    while IFS= read -rd '' dir; do 
      echo \"\\n============== working with $dir ===============\";
      grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=2e-4 --exclude-path \"/misc\" --exclude-path \"/git_revisions\" --exclude-path \"/MPI details\" \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"  \"${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/$dir/const.h5\"     &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/$dir/temp.xmf\"                                                          &&
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/$dir/timestep0000000000.h5\" rv      &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/$dir/timestep0000000000.h5\" th      &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_BINARY_DIR}/refdata_smg/$dir/timestep0000000000.xmf\"                              || exit 1; 
    done  
")

# Mark that these tests need the decompressed data
set_tests_properties(api_test_iles_diff api_test_smg_diff
                     PROPERTIES FIXTURES_REQUIRED "RefDataFixture")


# 3. Define the Cleanup Test (Delete the extracted folder)
add_test(NAME CleanupReferenceData 
         COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_BINARY_DIR}/refdata_iles ${CMAKE_CURRENT_BINARY_DIR}/refdata_smg)

# Mark this test as the cleanup step for the fixture
set_tests_properties(CleanupReferenceData PROPERTIES FIXTURES_CLEANUP "RefDataFixture")
