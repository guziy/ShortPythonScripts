from collections import OrderedDict
import os

__author__ = 'huziy'

DEFAULT_DATA_PATH = "/home/huziy/DATA/modis_surface_albedo"
from common import defaults
from pyhdf.SD import SD
import numpy as np

import mcd43c3_params

albedo_name = "Albedo_WSA_vis"

var_name_to_coef = {
    albedo_name: 1.0e-3
}

var_name_to_mask_value = {
    albedo_name: 32767
}

import biggus
from pyresample.geometry import SwathDefinition
from pyresample import image
from pyresample import kd_tree


class DataForDay(object):
    def __init__(self, var_name=albedo_name, day="2001.01.01"):
        self.shape = mcd43c3_params.data_shape
        self.dtype = np.dtype(np.float32)
        self.var_name = var_name
        self.path = os.path.join(DEFAULT_DATA_PATH, day)
        # select assuming there is only one hdf file in the folder
        for p in os.listdir(self.path):
            if p.endswith("hdf"):
                self.path = os.path.join(self.path, p)
                break

    def __getitem__(self, item):
        print self.path
        sd = SD(self.path)
        data = sd.select(self.var_name)[:]

        #data = np.flipud(data)
        data = np.ma.masked_where(data == var_name_to_mask_value[self.var_name],
                                  data)
        data = data.astype(self.dtype) * var_name_to_coef[self.var_name]

        return data[item]


def get_seasonal_mean_fields(season_name_to_months=None,
                             folder_path=DEFAULT_DATA_PATH,
                             start_year=2001, end_year=2010,
                             varname=albedo_name):
    if season_name_to_months is None:
        season_name_to_months = defaults.default_season_name_to_months

    all_folders = os.listdir(folder_path)
    all_years = [int(fname.split(".")[0]) for fname in all_folders]

    folder_names = [fname for fname, year in zip(all_folders, all_years) if start_year <= year <= end_year]
    months = [int(fname.split(".")[1]) for fname in folder_names]

    season_name_to_folder_list = OrderedDict([
        (sname, [fname for fname, m in zip(folder_names, months) if m in smonths])
        for sname, smonths in season_name_to_months.iteritems()
    ])

    result = OrderedDict()
    for sname, folders in season_name_to_folder_list.iteritems():
        arr_stack = biggus.ArrayStack(
            np.array([biggus.NumpyArrayAdapter(DataForDay(var_name=varname, day=fname)) for fname in folders])
        )
        print arr_stack.shape
        result[sname] = np.flipud(biggus.mean(arr_stack, axis=0).masked_array())

    return result


def interpolate_data_to(lons2d, lats2d, data, radius_of_influence_deg=0.5):
    ad = mcd43c3_params.pr_area_definition
    sd = SwathDefinition(lons=lons2d, lats=lats2d)

    wf = lambda r: 1.0 / r ** 2
    radius_of_influence = radius_of_influence_deg * np.pi / 180.0 * (6400.0 * 1000.0)  # convert to meters approximately

    nneighbours = int((radius_of_influence_deg / mcd43c3_params.dx) ** 2)
    print nneighbours

    return kd_tree.resample_custom(ad, data, sd, radius_of_influence, weight_funcs=wf, neighbours=nneighbours)


def main():
    from crcm5_canesm_mpi_comparisons.calculate_seasonal_means_and_significance \
        import get_basemap_and_coordinates_from_any_file_in

    basemap, lons2d, lats2d = get_basemap_and_coordinates_from_any_file_in()
    lons2d[lons2d > 180] -= 360



    data = get_seasonal_mean_fields(
        OrderedDict([("DJF", (1, 2, 12))]), end_year=2001
    )
    #print data["DJF"]
    data = data["DJF"]




    idata = interpolate_data_to(lons2d, lats2d, data)
    print idata.shape, idata.min(), idata.max(), idata.dtype

    import matplotlib.pyplot as plt

    plt.pcolormesh(data)




    plt.figure()
    plt.pcolormesh(idata)
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    main()