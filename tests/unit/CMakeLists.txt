add_executable(api_test api_test.cpp)

add_test(api_test_iles api_test ${CMAKE_BINARY_DIR})
add_test(api_test_smg  api_test ${CMAKE_BINARY_DIR} " --sgs=1 ")

add_test(api_test_iles_diff bash -c " 
  cd ${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"working with $dir\"                                                                                                                           
      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --exclude-path \"/git_revisions\" \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/const.h5\"                   &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/temp.xmf\"                                                         && 
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/timestep0000000000.h5\"                        && 
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_iles/$dir/timestep0000000000.xmf\"                             || exit 1; 
    done  
")

add_test(api_test_smg_diff bash -c " 
  cd ${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"working with $dir\"                                                                                                                           
      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --exclude-path \"/git_revisions\" \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/const.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/const.h5\"                   &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/temp.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/temp.xmf\"                                                         && 
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.h5\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/timestep0000000000.h5\"                        && 
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"${CMAKE_CURRENT_BINARY_DIR}/output/$dir/timestep0000000000.xmf\"     \"${CMAKE_CURRENT_SOURCE_DIR}/refdata_smg/$dir/timestep0000000000.xmf\"                             || exit 1; 
    done  
")
