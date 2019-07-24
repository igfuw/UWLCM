import numpy as np
import sys 
import os

def main():
    bin_type = sys.argv[1]

    wet_bins_str = gen_wet_bins()
    dry_bins_str = gen_dry_bins()

    if bin_type == "wet":
        print(wet_bins_str)
    elif bin_type == "dry":
        print(dry_bins_str)
    else:
        raise Exception('bin_type argument must be either "wet" or "dry"')

# generate 25 wet bins from 0.001 to 100 um
def gen_wet_bins():
    bin_str = ""
    left = 0.001e-6
    right = 100.0e-6
    nbins = 10
    moms = [0,1,3]
    moms_str = ','.join(str(x) for x in moms)
    str_fmt = "{:.1e}"
    left_edges = np.geomspace(left,right,num=nbins,endpoint=False,dtype=float)
    for i in np.arange(nbins):
        bin_str += str_fmt.format(left_edges[i])
        bin_str += ":" 
        if i < nbins-1:
            bin_str += str_fmt.format(left_edges[i+1])
            bin_str += "|"+moms_str+";"
        else:
            bin_str += str_fmt.format(right)
            bin_str += "|"+moms_str
    return bin_str 

# generate 40 dry bins from 0.001 to 10 um
def gen_dry_bins():
    bin_str = ""
    left = 0.001e-6
    right = 10.0e-6
    nbins = 10
    moms = [0,3]
    moms_str = ','.join(str(x) for x in moms)
    str_fmt = "{:.1e}"
    left_edges = np.geomspace(left,right,num=nbins,endpoint=False,dtype=float)
    for i in np.arange(nbins):
        bin_str += str_fmt.format(left_edges[i])
        bin_str += ":" 
        if i < nbins-1:
            bin_str += str_fmt.format(left_edges[i+1])
            bin_str += "|"+moms_str+";"
        else:
            bin_str += str_fmt.format(right)
            bin_str += "|"+moms_str
    return bin_str 

if __name__ == "__main__":
    main()