if(UWLCM_DISABLE)
  message(WARNING "Compilation of some of the options was disabled via the UWLCM_DISABLE flag. Unit tests are most likely to fail.") 
endif()

add_executable(api_test api_test.cpp)
target_compile_features(api_test PRIVATE cxx_std_11)

add_test(api_test_iles api_test ${CMAKE_BINARY_DIR})
add_test(api_test_smg  api_test ${CMAKE_BINARY_DIR} " --sgs=1 ")

add_test(api_test_iles_diff bash -c " 
  cd ${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"\\n============== working with $dir ===============\";
      grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=3e-5 --exclude-path \"/git_revisions\" --exclude-path \"/MPI details\" \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/const.h5\"                 &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/temp.xmf\"                                                         && 
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/timestep0000000000.h5\" rv     &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/timestep0000000000.h5\" th     &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/timestep0000000000.xmf\"                             || exit 1; 
    done  
")

add_test(api_test_smg_diff bash -c " 
  cd ${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"\\n============== working with $dir ===============\";
      grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=3e-5 --exclude-path \"/git_revisions\" --exclude-path \"/MPI details\" \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/const.h5\"                  &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/temp.xmf\"                                                          &&
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/timestep0000000000.h5\" rv      &&
      h5diff  -v2 --relative=1e-6 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/timestep0000000000.h5\" th      &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/timestep0000000000.xmf\"                              || exit 1; 
    done  
")
