#include "get_uwlcm_git_revision.hpp"
#include "../../../include/UWLCM/git_revision.h"

std::string get_uwlcm_git_revision()
{
  return UWLCM_GIT_REVISION;
}
