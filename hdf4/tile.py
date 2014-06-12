from datetime import datetime
from netCDF4 import Dataset
import os

__author__ = 'huziy'

from pyhdf.SD import SD
import numpy as np
import biggus

from scipy.ndimage import interpolation

# Tile image is here: http://nsidc.org/data/docs/daac/mod10_modis_snow/landgrid.html
#Documentation on tile bounds can be found here: #
#http://nsidc.org/data/docs/daac/mod10_modis_snow/sinusoidal_tile_coordinates.html

tile_hmax = 35  # goes from 0 to 35
tile_vmax = 17  # goes from 0 to 17


def get_area_properties_for(lower_left_indices=None, nx=None, ny=None):
    pass


class Tile(object):
    empty_tile_arr = None

    def __init__(self, h=0, v=0, data_folder="", date=None, varname="Snow_Albedo_Daily_Tile"):
        day_folder = "{0:%Y.%m.%d}".format(date)
        day_folder = os.path.join(data_folder, day_folder)

        #print h, v, "h{0:02d}v{1:02d}".format(h, v)

        fnames = [fname for fname in os.listdir(day_folder)
                  if fname.endswith("hdf") and "h{0:02d}v{1:02d}".format(h, v) in fname]
        self.empty_tile = False
        self.dtype = np.dtype("int32")
        self.varname = varname

        self.fill_value = -1

        if len(fnames) > 0:
            self.file_path = os.path.join(day_folder, fnames[0])
            ds = SD(self.file_path)
            hdf_var = ds.select(varname)

            self.shape = (hdf_var.dim(0).length(), hdf_var.dim(1).length())
            ds.end()
        else:
            self.shape = (2400, 2400)
            self.empty_tile = True


    def __getitem__(self, slices):
        if self.empty_tile:
            if self.empty_tile_arr is None:
                self.empty_tile_arr = -np.ones(self.shape, dtype=self.dtype)

            return self.empty_tile_arr[slices]
        ds = SD(self.file_path)
        data = ds.select(self.varname)[:]

        print data.dtype
        data = np.flipud(data)
        #data[(data < 0) | (data > 100)] = -1
        data = np.ma.masked_where((data < 0) | (data > 100), data)
        ds.end()
        return data[slices]


def calculate_seasonal_means_for_all_tiles(season_name=None, months=None, start_year=2001, end_year=2010,
                                           path="/b2_fs2/huziy/PythonProjects/DATA/modis/n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005"):
    for h in range(tile_hmax + 1):
        for v in range(tile_vmax + 1):
            #calculate_seasonal_mean_for_tile_without_biggus(h=h, v=v, season_name=season_name, path=path, months=months,
            #                                 start_year=start_year, end_year=end_year)

            if not (h == 14 and v == 3):
                continue

            calculate_seasonal_mean_for_tile(h=h, v=v, season_name=season_name, path=path, months=months,
                                             start_year=start_year, end_year=end_year)


def calculate_seasonal_mean_for_tile(h=0, v=0,
                                     path="/b2_fs2/huziy/PythonProjects/DATA/modis/n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005",
                                     start_year=2001, end_year=2010, season_name="DJF", months=None):
    dates = [datetime.strptime(folder_name, "%Y.%m.%d") for folder_name in os.listdir(path)]

    dates = [d for d in dates if (d.year >= start_year) and (d.year <= end_year) and (d.month in months)]

    #Tile(h=14,v=3,date=dates[0], data_folder=path)[:]
    
    #if True:
    #    return 


    arr_stack = biggus.ArrayStack(
        np.array([biggus.NumpyArrayAdapter(Tile(h=h, v=v, data_folder=path, date=d)) for d in dates]))

    the_mean = biggus.mean(arr_stack, axis=0).masked_array()

    folder_path = "/home/huziy/DATA/seasonal_modis_snow_albedo_biggus/{}".format(season_name)
    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)

    ds = Dataset(os.path.join(folder_path, "{}_h{}_v{}.nc".format(season_name, h, v)), mode="w")
    ds.createDimension("y", the_mean.shape[0])
    ds.createDimension("x", the_mean.shape[1])
    var_nc = ds.createVariable("I6", the_mean.dtype, ("y", "x"))
    var_nc.missing_value = -1
    the_mean[the_mean.mask] = var_nc.missing_value
    biggus.save([the_mean, ], [var_nc, ])

    ds.close()


def calculate_seasonal_mean_for_tile_without_biggus(h=0, v=0,
                                                    path="/b2_fs2/huziy/PythonProjects/DATA/modis/n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005",
                                                    start_year=2001, end_year=2010, months=None, season_name=""):
    dates = [datetime.strptime(folder_name, "%Y.%m.%d") for folder_name in os.listdir(path)]

    #select dates for the season
    dates = [d for d in dates if (d.year >= start_year) and (d.year <= end_year) and (d.month in months)]

    the_mean = Tile(date=dates.pop(), h=h, v=v, data_folder=path)[:]
    counter = 1.0
    import numexpr as ne

    for d in dates:
        print d, h, v
        current = Tile(date=d, h=h, v=v, data_folder=path)[:]
        the_mean = ne.evaluate("the_mean + current")
        counter += 1.0
    the_mean /= counter

    out_folder_path = "/home/huziy/DATA/seasonal_modis_snow_albedo/{}".format(season_name)

    if not os.path.isdir(out_folder_path):
        os.makedirs(out_folder_path)

    ds = Dataset(os.path.join(out_folder_path, "{}_h{}_v{}.nc".format(season_name, h, v)), mode="w")
    ds.createDimension("y", the_mean.shape[0])
    ds.createDimension("x", the_mean.shape[1])
    var_nc = ds.createVariable("I6", the_mean.dtype, ("y", "x"))
    var_nc.missing_value = -1
    var_nc[:] = the_mean
    #    the_mean[the_mean.mask] = var_nc.missing_value
    #    biggus.save([the_mean, ], [var_nc, ])

    ds.close()


def get_data_for_date(date=None,
                      path="/b2_fs2/huziy/PythonProjects/DATA/modis/n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005"):
    print "reading data for: {}".format(date)
    rows = []
    for j in range(tile_vmax + 1):
        row = [biggus.NumpyArrayAdapter(Tile(h=i, v=j, data_folder=path, date=date))
               for i in range(tile_hmax + 1)]
        row = biggus.LinearMosaic(row, axis=1)  # glueing along the second axes
        rows.append(row)

    big_arr = biggus.LinearMosaic(rows, axis=0)
    return big_arr


def interpolate_data_to():
    pass


def check_seasonal_mean():
    #calculate_seasonal_mean_for_tile(v=5, h=8, season_name="DJF", months=[1,2,12], start_year=2001, end_year=2010)
    calculate_seasonal_means_for_all_tiles(season_name="DJF", months=[1, 2, 12], start_year=2001, end_year=2010)
    # calculate_seasonal_means_for_all_tiles(season_name="MAM", months=[3, 4, 5], start_year=2001, end_year=2010)
    # calculate_seasonal_means_for_all_tiles(season_name="JJA", months=[6, 7, 8], start_year=2001, end_year=2010)
    # calculate_seasonal_means_for_all_tiles(season_name="SON", months=[9, 10, 11], start_year=2001, end_year=2010)


def check():
    path = "/b2_fs2/huziy/PythonProjects/DATA/modis/n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005"
    date = datetime(2000, 5, 24)
    rows = []
    for j in range(tile_vmax + 1):
        row = [biggus.NumpyArrayAdapter(Tile(h=i, v=j, data_folder=path, date=date)) for i in range(tile_hmax + 1)]
        row = biggus.LinearMosaic(row, axis=1)  # glueing along the second axes
        rows.append(row)

    big_arr = biggus.LinearMosaic(rows, axis=0)

    print big_arr.shape

    nx, ny = big_arr.shape

    import matplotlib.pyplot as plt

    plt.figure()
    subset = np.flipud(big_arr[::100, ::100].masked_array())

    print subset.min(), subset.max()
    #subset = np.ma.masked_where((subset > 100) | (subset < 0), subset)
    plt.pcolormesh(subset[:, :])
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    #check()
    check_seasonal_mean()
