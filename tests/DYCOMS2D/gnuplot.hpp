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
  gp << "set size square\n";
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
  gp << "set size square\n";
  gp << "set encoding utf8\n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";
  gp << "dx = 6400./" << n["x"]-1 << "\n"; 
  gp << "dy = 1500./" << n["z"]-1 << "\n"; 
  gp << "set xtics out scale .5 rotate by 60 ('0.0' 0, '.8' 800/dx, '1.6' 1600/dx, '2.4' 2400/dx, '3.2' 3200/dx, '4.0' 4000/dx, '4.8' 4800/dx, '5.6' 5600/dx, '6.4' 6400/dx)\n"; 
  gp << "set xlabel 'x [km]'\n";
  gp << "set ytics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dy, '0.6' 600/dy, '0.9' 900/dy, '1.2' 1200/dy, '1.5' 1500/dy)\n"; 
  gp << "set ylabel 'z [km]'\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
}

template <class data_t>
void plot(Gnuplot &gp, const data_t &data)
{
  blitz::Array<float, 2> tmp(data);

  gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
  gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
  gp << "splot '-' binary" << gp.binfmt(tmp.transpose(blitz::secondDim, blitz::firstDim)) << " scan=yx origin=(0,0,0) with image failsafe notitle\n";
  gp.sendBinary(tmp);
}
