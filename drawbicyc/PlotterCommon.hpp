#pragma once

#include <blitz/array.h>
#include <H5Cpp.h>
#include <map>


class PlotterCommon
{
  public:
  std::map<std::string, double> map;
  const string file;
  double CellVol, DomainSurf;

  protected:
  H5::H5File h5f;
  H5::DataSet h5d;
  H5::DataSpace h5s;

  void h5load(
    const string &file, 
    const string &dataset
  )
  {
    notice_macro("about to open file: " << file)
    h5f.openFile(file, H5F_ACC_RDONLY);

    notice_macro("about to read dataset: " << dataset)
    h5d = h5f.openDataSet(dataset);
    h5s = h5d.getSpace();
  }

  public:
  //ctor
  PlotterCommon(const string &file):
    file(file)
  {
    // init dt and outfreq
    {
      h5load(file + "/const.h5", "T");

      float dt;
      {
        auto attr = h5d.openAttribute("dt");
        attr.read(attr.getDataType(), &dt);
      }
      map["dt"] = dt;

      const hsize_t two = 2, zero = 0;
      float tmp[2];
      h5s.selectHyperslab( H5S_SELECT_SET, &two, &zero);
      h5d.read(tmp, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1, &two), h5s);
      map["outfreq"] = (tmp[1] - tmp[0]) / dt;
    }
  }
};

