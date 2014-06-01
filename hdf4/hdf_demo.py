__author__ = 'huziy'


import os
from rpn.rpn import RPN


from pyhdf.SD import SD
import pyhdf
import matplotlib.pyplot as plt

def main():

    path = "/b2_fs2/huziy/PythonProjects/DATA/modis_test/MOD10A1.A2001145.h23v01.005.2007024165028.hdf"

    vname = "Snow_Albedo_Daily_Tile"
    ds = SD(path)
    print ds.datasets()


    albedo = ds.select(vname)
    assert isinstance(albedo, pyhdf.SD.SDS)
    print albedo.info()
    print albedo.getrange()
    print albedo.iscoordvar()
    print albedo.dimensions()
    print albedo.attributes()
    print albedo.ref()
    print albedo.dim(0)

    #var_data = ds.select(varname)[0, :, :]









if __name__ == "__main__":
    main()