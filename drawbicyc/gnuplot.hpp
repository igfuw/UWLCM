#pragma once

#include <blitz/array.h>
#include <gnuplot-iostream.h>
#include <map>


void init_prof(
  Gnuplot &gp, 
  const std::string &file, 
  const int &ny, const int &nx
)
{
  boost::filesystem::create_directories(
    boost::filesystem::path(file).parent_path()
  );
  gp << "set term svg dynamic enhanced fsize 16 size " << nx * 500 << "," << ny * 500 << "\n";
//  gp << "set size square\n";
//  gp << "set size ratio 0.3\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
  gp << "set yrange[0:1.2]\n";
  gp << "set border lw 1.5\n";
  gp << "set linetype 1 lw 2\n";
  gp << "set linetype 2 lw 2\n";
  gp << "set linetype 3 lw 2\n";
  gp << "set linetype 4 lw 2 lc rgb 'violet'\n";
  gp << "set linetype 5 lw 2\n";
  gp << "set linetype 6 lw 2\n";
}

void init(
  Gnuplot &gp, 
  const std::string &file, 
  const int &ny, const int &nx, 
  std::map<std::string, double> n
)
{
  boost::filesystem::create_directories(
    boost::filesystem::path(file).parent_path()
  );

  gp << "set term svg dynamic enhanced fsize 13 size " << nx * 500 << "," << ny * 500 << "\n";
//  gp << "set size square\n";
  gp << "set size ratio 0.66666666 \n";
  gp << "set encoding utf8\n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";
  gp << "dx = 3600./" << n["x"]-1 << "\n";   // TODO: case-specific
  gp << "dy = 2400./" << n["z"]-1 << "\n";   // TODO: case-specific
  gp << "set xtics out scale .5 rotate by 60 ('0.0' 0, '0.6' 600/dx, '1.2' 1200/dx, '1.8' 1800/dx, '2.4' 2400/dx, '3.0' 3000/dx, '3.6' 3600/dx)\n"; 
  gp << "set xlabel 'x [km]'\n";
  gp << "set ytics out scale .5 rotate by 60 ('0.0' 0, '0.6' 600/dy, '1.2' 1200/dy, '1.8' 1800/dy, '2.4' 2400/dy)\n"; 
  gp << "set ylabel 'z [km]'\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
}
