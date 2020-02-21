# CMake generated Testfile for 
# Source directory: /home/pzmij/biblioteki/UWLCM/tests/unit
# Build directory: /home/pzmij/biblioteki/UWLCM/build/tests/unit
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(api_test_iles "api_test" "/home/pzmij/biblioteki/UWLCM/build")
add_test(api_test_smg "api_test" "/home/pzmij/biblioteki/UWLCM/build" " --sgs=1 ")
add_test(api_test_iles_diff "bash" "-c" " 
  cd /home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"\n============== working with $dir ===============\";
      grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=3e-5 --exclude-path \"/git_revisions\" \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/const.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/$dir/const.h5\"                 &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/temp.xmf\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/$dir/temp.xmf\"                                                         && 
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/$dir/timestep0000000000.h5\" rv     &&
      h5diff  -v2 --relative=1e-6 \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/$dir/timestep0000000000.h5\" th     &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.xmf\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_iles/$dir/timestep0000000000.xmf\"                             || exit 1; 
    done  
")
add_test(api_test_smg_diff "bash" "-c" " 
  cd /home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/ && find ./* -maxdepth 1 -type d -print0 |            
    while IFS= read -rd '' dir; do 
      echo \"\n============== working with $dir ===============\";
      grep $(echo \"$dir\" | sed 's/..//') hash_dict.txt

      echo   'comparing const.h5'                                                                                                                          
      h5diff -v2 --relative=3e-5 --exclude-path \"/git_revisions\" \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/const.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/$dir/const.h5\"                  &&
      echo   'comparing temp.xmf'                                                                                                                                                      &&
      diff \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/temp.xmf\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/$dir/temp.xmf\"                                                          &&
      echo   'comparing timestep0000000000.h5'                                                                                                                                         &&
      h5diff  -v2 --relative=1e-6 \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/$dir/timestep0000000000.h5\" rv      &&
      h5diff  -v2 --relative=1e-6 \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.h5\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/$dir/timestep0000000000.h5\" th      &&
      echo   'comparing timestep0000000000.xmf'                                                                                                                                        &&
      diff \"/home/pzmij/biblioteki/UWLCM/build/tests/unit/output/$dir/timestep0000000000.xmf\"     \"/home/pzmij/biblioteki/UWLCM/tests/unit/refdata_smg/$dir/timestep0000000000.xmf\"                              || exit 1; 
    done  
")
