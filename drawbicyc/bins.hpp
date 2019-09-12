#pragma once

// bin sizes for calc and plot
/*
vector<quantity<si::length> > bins_dry()
{
  vector<quantity<si::length> > ret;
  // dry radius bins: .001 ... .01 ... 10 (40 bins in total)
  for (int i = 0; i < 40; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .1) * si::metres);
  return ret;
}
*/

vector<quantity<si::length> > bins_wet()
{
  vector<quantity<si::length> > ret;
  for (int i = 0; i < 16; ++i)
    ret.push_back(2e-6 * i * si::metres); 
  return ret;
}

// focus plot locations
const int yslice_idx=75;
const int z_idx=100;
/*
std::pair<
  std::set<std::pair<int,int> >,
  std::set<std::pair<int,int> >
> focus_2d = {
  {   // left column
    {16+ox, 20+oz}, 
    {16+ox, 37+oz}, 
    {16+ox, 39+oz}, 
    {16+ox, 44+oz},
    {16+ox, 65+oz} 
  },{ // right column
    {58+ox, 20+oz},
    {58+ox, 37+oz}, 
    {58+ox, 39+oz}, 
    {58+ox, 44+oz},
    {58+ox, 65+oz}
  }
};
*/

std::set<std::array<int, 3>> focus_3d = {
  {60, yslice_idx, z_idx }, 
  {68, yslice_idx, z_idx }, 
  {76, yslice_idx, z_idx },
  {84, yslice_idx, z_idx },
  {92, yslice_idx, z_idx } 
};

const int box_size=2; // range of cells in each direction in the box is [x-box_size, x+box_size]
