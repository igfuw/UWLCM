#pragma once

#include <iomanip> // setprecision
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
  gp << "set term svg dynamic enhanced size " << nx * 500 << "," << ny * 500 << " font ',16'\n";
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
  gp << "set nokey\n";
}

void init(
  Gnuplot &gp, 
  const std::string &file, 
  const int &ny, const int &nx, // number of multiplots in y and x 
  std::map<std::string, double> n,
  double size_scale = 1,
  double ratio = 0 // default used to be 0.666666666
)
{
  boost::filesystem::create_directories(
    boost::filesystem::path(file).parent_path()
  );

  if(ratio == 0)
    ratio = n["dz"] / n["dx"];

  const int xtics = 7;
  const int ytics = xtics * ratio + 0.5;

  gp << "set term svg dynamic enhanced size " << nx * size_scale * 500 << "," << ny * size_scale * 500 << " font ',13'\n";
//  gp << "set size square\n";
  gp << "set size ratio "<< ratio <<" \n";
  gp << "set encoding utf8\n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";
  gp << "dx = " << n["dx"] << "\n"; 
  gp << "dy = " << n["dy"] << "\n"; 
  gp << "dz = " << n["dz"] << "\n"; 

  gp << "set xtics out scale .5 rotate by 60 (";
  for(int i=0; i<xtics; ++i)
  {
    double label = double(i) / (xtics-1)* n["x"]  * n["dx"] / 1e3; 
    gp << "'" << std::setprecision(1) << label << "' " << double(i * n["x"] / (xtics-1));
    if(i < xtics-1)
      gp << ", ";
  } 
  gp << ")\n"; 

  gp << "set ytics out scale .5 rotate by 60 (";
  for(int i=0; i<ytics; ++i)
  {
    double label = double(i) / (ytics-1)* n["z"]  * n["dz"] / 1e3; 
    gp << "'" << std::setprecision(1) << label << "' " << double(i * n["z"] / (ytics-1));
    if(i < ytics-1)
      gp << ", ";
  } 
  gp << ")\n"; 

  gp << "set xlabel 'x [km]'\n";
  gp << "set ylabel 'z [km]'\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
}
