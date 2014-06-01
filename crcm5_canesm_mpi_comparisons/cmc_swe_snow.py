from datetime import datetime
from netCDF4 import Dataset, date2num
import os
import re
import collections
import numpy as np
from scipy.spatial.ckdtree import cKDTree
from geo import lat_lon


__author__ = 'huziy'


class CmcSweSnowManager(object):
    def __init__(self, path_snow_depth="", path_swe=""):
        """
        First dimension of the data corresponds to longitude and the second to latitude
        :param path_snow_depth:
        :param path_swe:
        """
        self.path_snow_depth, self.path_swe = get_default_snow_depth_and_swe_paths()

        if len(path_snow_depth) > 0:
            self.path_snow_depth = path_snow_depth

        if len(path_swe) > 0:
            self.path_swe = path_swe

        self.lons2d = None
        self.lats2d = None
        self._read_coordinates()

        self.land_sea_mask = None
        self._read_land_sea_mask()

        self.swe_varname = "swe"
        self.sd_varname = "sd"

        self.kdtree = None

    def interpolate_data_to(self, target_lons=None, target_lats=None, data_field=None, nneighbours=1):
        if self.kdtree is None:
            xs, ys, zs = lat_lon.lon_lat_to_cartesian(self.lons2d.flatten(), self.lats2d.flatten())
            self.kdtree = cKDTree(zip(xs, ys, zs))

        xt, yt, zt = lat_lon.lon_lat_to_cartesian(target_lons.flatten(), target_lats.flatten())
        dists, inds = self.kdtree.query(zip(xt, yt, zt), k=nneighbours)


        data_confident = np.ma.masked_where(self.land_sea_mask, data_field)

        if nneighbours == 1:
            return data_field.flatten()[inds].reshape(target_lons.shape)

        w = 1.0 / dists ** 2
        nearest_data = data_confident.flatten()[inds]
        wx = nearest_data * w



        result = wx.sum(axis=1) / w.sum(axis=1)
        result = np.ma.masked_where(nearest_data.mask.astype(int).sum(axis=1) / float(w.shape[1]) > 0.5, result)
        return result.reshape(target_lons.shape)




    def convert_data_to_netcdf(self, outfile_swe="/home/huziy/DATA/CMC/cmc_swe.nc",
                               outfile_sd="/home/huziy/DATA/CMC/cmc_sd.nc"):
        missing_value = -999.0

        ds = Dataset(outfile_swe, "w", format="NETCDF3_CLASSIC")

        ds.createDimension("time")
        ny, nx = self.lons2d.shape
        ds.createDimension("lat", ny)
        ds.createDimension("lon", nx)

        lon_var = ds.createVariable("longitude", "f4", ("lat", "lon"))
        lat_var = ds.createVariable("latitude", "f4", ("lat", "lon"))
        lon_var[:] = self.lons2d
        lat_var[:] = self.lats2d

        #define time for the swe variable
        time_swe_var = ds.createVariable("time", "i4", ("time", ))
        time_swe_var.units = "days since 1997-01-01 00:00:00"


        #create swe variable
        swe_var = ds.createVariable("swe", "f4", ("time", "lat", "lon"))
        swe_var.units = "mm"
        swe_var.description = "Monthly mean snow water equivalent"
        swe_var.missing_value = missing_value


        #read and save swe data
        all_dates = []
        year_to_month_to_swe = self._read_all_monthly_swe_data()

        record = 0
        for year in sorted(year_to_month_to_swe):
            month_to_swe = year_to_month_to_swe[year]
            for month in sorted(month_to_swe):
                arr = np.asarray(month_to_swe[month])
                arr[self.land_sea_mask] = missing_value
                swe_var[record, :, :] = arr
                record += 1
                all_dates.append(datetime(year, month, 15))

        time_swe_var[:] = date2num(all_dates, time_swe_var.units)
        ds.close()


        #Save snow depth to a different file
        ds = Dataset(outfile_sd, "w", format="NETCDF3_CLASSIC")
        ds.createDimension("time")
        ny, nx = self.lons2d.shape
        ds.createDimension("lat", ny)
        ds.createDimension("lon", nx)

        lon_var = ds.createVariable("longitude", "f4", ("lat", "lon"))
        lat_var = ds.createVariable("latitude", "f4", ("lat", "lon"))
        lon_var[:] = self.lons2d
        lat_var[:] = self.lats2d


        #create snow depth variable
        sd_var = ds.createVariable("snow_depth", "f4", ("time", "lat", "lon"))
        sd_var.units = "cm"
        sd_var.description = "Monthly mean snow depth"
        sd_var.missing_value = missing_value


        #define time for snow depth variable
        time_sd_var = ds.createVariable("time", "i4", ("time", ))
        time_sd_var.units = "days since 1997-01-01 00:00:00"


        #read and save snow depth data
        all_dates = []
        year_to_month_to_sd = self._read_all_monthly_snow_depth_data()

        record = 0
        for year in sorted(year_to_month_to_sd):
            month_to_sd = year_to_month_to_sd[year]
            for month in sorted(month_to_sd):
                arr = np.asarray(month_to_sd[month])
                arr[self.land_sea_mask] = missing_value
                sd_var[record, :, :] = arr
                record += 1
                all_dates.append(datetime(year, month, 15))

        time_sd_var[:] = date2num(all_dates, time_sd_var.units)

        ds.close()


    def _read_coordinates(self):
        fname = "cmc_analysis_ps_lat_long.txt"
        fpath = os.path.join(self.path_snow_depth, "..", fname)

        read_data = False
        with open(fpath) as f:
            for line in f:
                if "ni" in line and "nj" in line:
                    ni, nj = [int(g) for g in re.findall(r"\d+", line)]
                    self.lons2d = np.zeros((ni, nj))
                    self.lats2d = np.zeros((ni, nj))

                if read_data:
                    fields = [field.strip() for field in line.split()]
                    i, j = [int(g) - 1 for g in fields[:2]]
                    lat, lon = [float(g) for g in fields[2:]]
                    self.lons2d[i, j] = lon
                    self.lats2d[i, j] = lat

                if not read_data:
                    read_data = line.strip().startswith("----------")

    def get_land_sea_mask_interpolated_to(self, lons2d_target, lats2d_target):

        if self.kdtree is None:
            pass

        raise NotImplementedError()


    def _read_snow_depth_data_for_year(self, year):
        #returns a dictionary month => field

        fnames = [name for name in os.listdir(self.path_snow_depth) if str(year) in name and name.endswith(".txt")]
        assert len(fnames) == 1
        fname = fnames[0]

        fpath = os.path.join(self.path_snow_depth, fname)

        month_to_field = {}
        month = None
        with open(fpath) as f:
            for line in f:
                if line == "":
                    continue

                fields = line.split()
                if len(fields) == 2:
                    month = int(line.split()[-1].strip())
                    if month not in month_to_field:
                        month_to_field[month] = []
                else:
                    assert month is not None
                    month_to_field[month].append([float(field.strip()) for field in fields])

        for month in month_to_field:
            month_to_field[month] = np.asarray(month_to_field[month])
        return month_to_field

    def _read_all_monthly_snow_depth_data(self):
        fnames = [name for name in os.listdir(self.path_snow_depth) if name.endswith(".txt")]

        year_to_month_to_field = {}
        for fname in fnames:

            fpath = os.path.join(self.path_snow_depth, fname)

            month_to_field = {}
            month = None
            with open(fpath) as f:
                for line in f:
                    if line == "":
                        continue

                    fields = line.split()
                    if len(fields) == 2:
                        month = int(line.split()[-1].strip())
                        if month not in month_to_field:
                            month_to_field[month] = []
                    else:
                        assert month is not None
                        month_to_field[month].append([float(field.strip()) for field in fields])

            for month in month_to_field:
                month_to_field[month] = np.asarray(month_to_field[month])

            year = int(os.path.splitext(fname)[0].split("_")[-1])
            year_to_month_to_field[year] = month_to_field
        return year_to_month_to_field


    def _read_all_monthly_swe_data(self):
        fname = [name for name in os.listdir(self.path_swe) if name.endswith(".txt")][0]  # there is only one file

        fpath = os.path.join(self.path_swe, fname)
        with open(fpath) as f:
            current_month = None
            current_year = None
            year_to_month_to_field = {}
            for line in f:
                if line.strip() == "":
                    continue

                fields = line.split()

                if len(fields) == 3:  # New month
                    the_year, the_month, _ = [int(the_field) for the_field in fields]
                    current_month = the_month
                    if current_year != the_year:
                        year_to_month_to_field[the_year] = {}
                        current_year = the_year

                    year_to_month_to_field[current_year][current_month] = []
                else:

                    year_to_month_to_field[current_year][current_month].append(
                        [float(the_field) for the_field in fields]
                    )

        return year_to_month_to_field


    def get_seasonal_mean_clim(self, start_year=1999, end_year=2012,
                               season_name_to_months=None, var_name="sd"):
        """
        return {season: field_snow_depth}, for var_name specified
        varname can be either "swe" or "sd" -- snow water equivalent or snow depth respectively

        :param start_year:
        :param end_year:
        :param season_name_to_months:
        :raise ValueError:
        """
        min_start_year = 1998
        if start_year == min_start_year:
            print "Warning it is better not to use {} since the data for this year starts in April".format(start_year)

        if start_year < min_start_year:
            raise ValueError("No data available before {}".format(min_start_year))

        if var_name == "sd":
            year_to_month_to_data = self._read_all_monthly_snow_depth_data()
        else:
            year_to_month_to_data = self._read_all_monthly_swe_data()

        result = collections.OrderedDict()

        for season, month_list in season_name_to_months.iteritems():
            result[season] = np.mean([
                year_to_month_to_data[y][m] for y in range(start_year, end_year + 1) for m in month_list
            ], axis=0)

        return result



    def _read_land_sea_mask(self):
        fname = "cmc_analysis_lsmask_binary_nogl.txt"
        fpath = os.path.join(self.path_snow_depth, "..", fname)
        bits = []
        with open(fpath) as f:
            for line in f:
                if line.strip() == "":
                    continue
                bits.append([int(b) for b in re.findall(r"\d", line)])

        self.land_sea_mask = np.asarray(bits) < 0.5


def get_default_snow_depth_and_swe_paths():
    path_snow_depth = "/b2_fs2/huziy/PythonProjects/DATA/CMC/nsidc0447_CMC_snow_depth_v01/Monthly_Snow_Depth_Files"
    path_swe = "/b2_fs2/huziy/PythonProjects/DATA/CMC/nsidc0447_CMC_snow_depth_v01/Monthly_SWE_Estimates"
    return path_snow_depth, path_swe


def check():
    cmc = CmcSweSnowManager()
    correct_dims = (706, 706)
    assert cmc.lons2d.shape == correct_dims, "Dimensions read from the file is wrong: {}x{}, should be {}x{} ".format(
        *(cmc.lons2d.shape + correct_dims)
    )

    assert cmc.land_sea_mask.shape == correct_dims
    assert cmc.land_sea_mask.min() == 0
    assert cmc.land_sea_mask.max() == 1

    cmc.convert_data_to_netcdf()

    #cmc.get_seasonal_mean_clim()
    #import matplotlib.pyplot as plt

    #plt.figure()
    #plt.pcolormesh(cmc.land_sea_mask)
    #plt.show()



    print "You are awesome!!"
    print "All tests passed!!"


if __name__ == "__main__":
    check()