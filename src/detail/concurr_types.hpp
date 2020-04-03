#pragma once

#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/bcond/rigid_2d.hpp>
#include <libmpdata++/bcond/open_2d.hpp>
#include <libmpdata++/bcond/gndsky_2d.hpp>

#include <libmpdata++/bcond/cyclic_3d.hpp>
#include <libmpdata++/bcond/rigid_3d.hpp>
#include <libmpdata++/bcond/open_3d.hpp>
#include <libmpdata++/bcond/gndsky_3d.hpp>

#include <libmpdata++/concurr/openmp.hpp>

template <class solver_t, int n_dims>
struct concurr_openmp_rigid;

// 2D
template <class solver_t>
struct concurr_openmp_rigid<solver_t, 2>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::rigid,  libmpdataxx::bcond::rigid,
    libmpdataxx::bcond::rigid,  libmpdataxx::bcond::rigid 
  >;
};

// 3D
template <class solver_t>
struct concurr_openmp_rigid<solver_t, 3>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::rigid,  libmpdataxx::bcond::rigid,
    libmpdataxx::bcond::rigid,  libmpdataxx::bcond::rigid,
    libmpdataxx::bcond::rigid,  libmpdataxx::bcond::rigid 
  >;
};

template <class solver_t, int n_dims>
struct concurr_openmp_cyclic_gndsky;

// 2D
template <class solver_t>
struct concurr_openmp_cyclic_gndsky<solver_t, 2>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::gndsky, libmpdataxx::bcond::gndsky
  >;
};

// 3D
template <class solver_t>
struct concurr_openmp_cyclic_gndsky<solver_t, 3>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::gndsky, libmpdataxx::bcond::gndsky
  >;
};

template <class solver_t, int n_dims>
struct concurr_openmp_cyclic;

// 2D
template <class solver_t>
struct concurr_openmp_cyclic<solver_t, 2>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic
  >;
};

// 3D
template <class solver_t>
struct concurr_openmp_cyclic<solver_t, 3>
{
  using type = libmpdataxx::concurr::openmp<
    solver_t, 
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic,
    libmpdataxx::bcond::cyclic, libmpdataxx::bcond::cyclic
  >;
};
