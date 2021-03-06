import calendar
from datetime import timedelta, datetime
import itertools
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import pandas
from scipy.spatial.ckdtree import cKDTree
from scipy.spatial.kdtree import KDTree

from geo import lat_lon

__author__ = 'huziy'

import numpy as np
from netCDF4 import Dataset, num2date, date2num
import matplotlib.pyplot as plt
from collections import OrderedDict


class CRUDataManager:
    def __init__(self, path="data/cru_data/CRUTS3.1/cru_ts_3_10.1901.2009.tmp.dat.nc", var_name="tmp", lazy=False):
        self.times = None
        self.var_data = None

        self.times_var = None
        self.kdtree = None
        self.times_num = None
        self.lons2d, self.lats2d = None, None

        self.lazy = lazy
        ds = Dataset(path)
        self.var_name = var_name
        self._init_fields(ds)
        self.nc_dataset = ds

        self.nc_vars = ds.variables

        pass

    def _init_fields(self, nc_dataset):
        nc_vars = nc_dataset.variables
        lons = nc_vars["lon"][:]
        lats = nc_vars["lat"][:]

        lats2d, lons2d = np.meshgrid(lats, lons)

        self.lons2d, self.lats2d = lons2d, lats2d

        self.times_var = nc_vars["time"]
        self.times_num = nc_vars["time"][:]
        self.times = num2date(self.times_num, self.times_var.units, self.times_var.calendar)
        if not self.lazy:
            self.var_data = np.transpose(nc_vars[self.var_name][:], axes=[0, 2, 1])

        x_in, y_in, z_in = lat_lon.lon_lat_to_cartesian(self.lons2d.flatten(), self.lats2d.flatten())
        self.kdtree = cKDTree(zip(x_in, y_in, z_in))


    def get_seasonal_means(self, season_name_to_months=None, start_year = None, end_year = None):
        if season_name_to_months is None:
            season_name_to_months = OrderedDict([
                ("Winter", (1, 2, 12)),
                ("Spring", range(3, 6)),
                ("Summer", range(6, 9)),
                ("Fall", range(9, 12))])



        season_name_to_coef = {}
        for sname, months in season_name_to_months.iteritems():
            season_name_to_coef[sname] = 1

            if self.var_name.lower() == "pre":
                days = sum([calendar.monthrange(y, m)[1] for m in months for y in range(start_year, end_year + 1)])
                season_name_to_coef[sname] = 1.0 / float(days)


        month_to_season = {}
        for sname, mlist in season_name_to_months.iteritems():
            for m in mlist:
                month_to_season[m] = sname

        if self.var_data is None:
            self.var_data = np.transpose(self.nc_dataset.variables[self.var_name][:], axes=[0, 2, 1])

        nt, nx, ny = self.var_data.shape
        panel = pandas.Panel(data=self.var_data, items=self.times, major_axis=range(nx), minor_axis=range(ny))
        panel = panel.select(lambda d: start_year <= d.year <= end_year)

        if self.var_name == "pre":
            panel_seasonal = panel.groupby(lambda d: month_to_season[d.month], axis = "items").sum()
        else:
            panel_seasonal = panel.groupby(lambda d: month_to_season[d.month], axis = "items").mean()

        season_to_mean = OrderedDict()
        for sname, _ in season_name_to_months.iteritems():
            season_to_mean[sname] = panel_seasonal[sname].values * season_name_to_coef[sname]

        return season_to_mean


    def get_mean(self, start_year, end_year, months=None):
        """
        returns the mean for the period [start_year, end_year], over the months
        :type months: list
        months = list of month numbers over which the averaging is done
        """

        if months is None:
            months = list(range(1, 13))

        start_date = datetime(start_year, 1, 1)
        end_date = datetime(end_year + 1, 1, 1)

        start_date_num = date2num(start_date, self.times_var.units)
        end_date_num = date2num(end_date, self.times_var.units)

        sel_query = (self.times_num >= start_date_num) & (self.times_num < end_date_num)
        sel_dates = self.times_num[sel_query]
        sel_data = np.transpose(self.nc_vars[self.var_name][sel_query, :, :], axes=[0, 2, 1])

        sel_dates = num2date(sel_dates, self.times_var.units)

        ind_vector = np.where(map(lambda x: (x.month in months), sel_dates))[0]
        return np.mean(sel_data[ind_vector, :, :], axis=0)


    def get_daily_climatology(self, start_year, end_year, stamp_year=2001):
        """
        returns a numpy array of shape (365, nx, ny) with daily climatological means
        """
        day = timedelta(days=1)
        the_date = datetime(stamp_year, 1, 1)
        stamp_days = [the_date + i * day for i in xrange(365)]
        result = []

        nt, nx, ny = self.var_data.shape
        data_panel = pandas.Panel(data=self.var_data, items=self.times, major_axis=range(nx), minor_axis=range(ny))
        data_panel = data_panel.select(
            lambda d: (start_year <= d.year <= end_year) and not (d.day == 29 and d.month == 2))

        data_panel = data_panel.groupby(lambda d: datetime(stamp_year, d.month, d.day), axis="items").mean()
        assert isinstance(data_panel, pandas.Panel)
        data_panel = data_panel.sort_index()
        print data_panel.values.shape
        print data_panel.items


        #for the_date in stamp_days:
        #    bool_vector = np.array(map(lambda x: (x.day == the_date.day) and
        #                                         (x.month == the_date.month) and
        #                                         (x.year <= end_year) and (x.year >= start_year), self.times))
        #    result.append(np.mean(self.var_data[bool_vector, :, :], axis=0))
        #return np.array(result)
        return data_panel.values


    def interpolate_daily_climatology_to(self, clim_data, lons2d_target=None, lats2d_target=None):
        #expects clim_data to have the following shape (365, nx, ny)
        #        lons2d_target: (nx, ny)
        #        lats2d_target: (nx, ny)


        x, y, z = lat_lon.lon_lat_to_cartesian(lons2d_target.flatten(), lats2d_target.flatten())

        nt = clim_data.shape[0]
        data_help = np.reshape(clim_data, (nt, -1))

        dists, inds = self.kdtree.query(zip(x, y, z))

        return data_help[:, inds].reshape((nt,) + lons2d_target.shape)

        pass


    def get_thawing_index_from_climatology(self, daily_temps_clim, t0=0.0):

        nt, nx, ny = daily_temps_clim.shape
        result = np.zeros((nx, ny))

        for t in xrange(nt):
            tfield = daily_temps_clim[t, :, :]
            result += tfield * np.array(tfield >= t0).astype(int)
        return result

        pass


    def create_monthly_means_file(self, start_year, end_year):
        fname = "{0}_monthly_means.nc".format(self.var_name)
        year_range = range(start_year, end_year + 1)
        dsm = Dataset(fname, "w", format="NETCDF3_CLASSIC")
        dsm.createDimension('year', len(year_range))
        dsm.createDimension("month", 12)
        dsm.createDimension('lon', self.lons2d.shape[0])
        dsm.createDimension('lat', self.lons2d.shape[1])

        lonVariable = dsm.createVariable('longitude', 'f4', ('lon', 'lat'))
        latVariable = dsm.createVariable('latitude', 'f4', ('lon', 'lat'))
        yearVariable = dsm.createVariable("year", "i4", ("year",))

        variable = dsm.createVariable(self.var_name, "f4", ('year', "month", 'lon', 'lat'))
        for i, the_year in enumerate(year_range):
            print the_year
            for j, the_month in enumerate(xrange(1, 13)):
                variable[i, j, :, :] = self.get_mean(the_year, the_year, months=[the_month])

        lonVariable[:] = self.lons2d
        latVariable[:] = self.lats2d
        yearVariable[:] = np.array(year_range)
        dsm.close()

        pass


    def _interp_and_sum(self, data1d, mults_1d, x, y, z, nneighbors=1):
        data_interp = self.interpolate_data_to_cartesian(data1d, x, y, z, nneighbours=nneighbors)
        return np.sum(mults_1d * data_interp)




    def get_mean_upstream_timeseries_monthly(self, model_point, data_manager):
        """
        get mean swe upstream of the model_point

        year range for selection is in model_point.continuous_data_years() ..
        """


        #create the mask of points over which the averaging is going to be done
        lons_targ = data_manager.lons2D[model_point.flow_in_mask == 1]
        lats_targ = data_manager.lats2D[model_point.flow_in_mask == 1]

        xt, yt, zt = lat_lon.lon_lat_to_cartesian(lons_targ, lats_targ)

        nxs, nys = self.lons2d.shape
        i_source, j_source = range(nxs), range(nys)

        j_source, i_source = np.meshgrid(j_source, i_source)

        i_source = i_source.flatten()
        j_source = j_source.flatten()

        dists, inds = self.kdtree.query(zip(xt, yt, zt), k=1)
        ixsel = i_source[inds]
        jysel = j_source[inds]

        print "Calculating spatial mean"
        #calculate spatial mean
        #calculate spatial mean
        if self.lazy:
            theVar = self.nc_vars[self.var_name]

            data_series = []
            for i, j in zip(ixsel, jysel):
                data_series.append(theVar[:, j, i])

            data_series = np.mean(data_series, axis=0)
        else:
            data_series = np.mean(self.var_data[:, ixsel, jysel], axis=1)

        print "Finished calculating spatial mean"

        #calculate daily climatology
        df = pandas.DataFrame(data=data_series, index=self.times, columns=["values"])

        df["year"] = df.index.map(lambda d: d.year)

        df = df[df["year"].isin(model_point.continuous_data_years)]
        monthly_clim = df.groupby(by=lambda d: d.month).mean()

        month_dates = [datetime(1985, m, 15) for m in range(1, 13)]
        vals = [monthly_clim.ix[d.month, "values"] for d in month_dates]

        return pandas.TimeSeries(data=vals, index=month_dates)


    def get_mean_upstream_timeseries_daily(self, model_point, dm, stamp_dates=None):
        """
        get mean swe upstream of the model_point
        """


        #create the mask of points over which the averaging is going to be done
        lons_targ = dm.lons2D[model_point.flow_in_mask == 1]
        lats_targ = dm.lats2D[model_point.flow_in_mask == 1]

        xt, yt, zt = lat_lon.lon_lat_to_cartesian(lons_targ, lats_targ)

        nxs, nys = self.lons2d.shape
        i_source, j_source = range(nxs), range(nys)

        j_source, i_source = np.meshgrid(j_source, i_source)

        i_source = i_source.flatten()
        j_source = j_source.flatten()

        dists, inds = self.kdtree.query(zip(xt, yt, zt), k=1)
        ixsel = i_source[inds]
        jysel = j_source[inds]

        df_empty = pandas.DataFrame(index=self.times)
        df_empty["year"] = df_empty.index.map(lambda d: d.year)

        #calculate spatial mean
        #calculate spatial mean
        sel_date_indices = np.where(df_empty["year"].isin(model_point.continuous_data_years))[0]
        if self.lazy:
            theVar = self.nc_vars[self.var_name]
            data_series = np.mean([theVar[sel_date_indices, j, i] for i, j in zip(ixsel, jysel)], axis=0)
        else:
            data_series = np.mean(self.var_data[:, ixsel, jysel], axis=1)


        #calculate daily climatology
        df = pandas.DataFrame(data=data_series, index=self.times, columns=["values"])

        df["year"] = df.index.map(lambda d: d.year)
        df = df[df["year"].isin(model_point.continuous_data_years)]
        daily_clim = df.groupby(by=lambda d: (d.month, d.day)).mean()

        vals = [daily_clim.ix[(d.month, d.day), "values"] for d in stamp_dates]
        return pandas.TimeSeries(data=vals, index=stamp_dates)

    def interpolate_data_to_cartesian(self, data_in_flat, x, y, z, nneighbours=4):
        """
        len(data_in_flat) , len(x) == len(y) == len(z) == len(data_out_flat) - all 1D
        """
        print "start query"
        dst, ind = self.kdtree.query(zip(x, y, z), k=nneighbours)
        print "end query"

        inverse_square = 1.0 / dst ** 2
        if len(dst.shape) > 1:
            norm = np.sum(inverse_square, axis=1)
            norm = np.array([norm] * dst.shape[1]).transpose()
            coefs = inverse_square / norm

            data_out_flat = np.sum(coefs * data_in_flat[ind], axis=1)
        elif len(dst.shape) == 1:
            data_out_flat = data_in_flat[ind]
        else:
            raise Exception("Could not find neighbor points")
        return data_out_flat


    def interpolate_data_to(self, data_in, lons2d, lats2d, nneighbours=4):
        """
        Interpolates data_in to the grid defined by (lons2d, lats2d)
        assuming that the data_in field is on the initial CRU grid

        interpolate using 4 nearest neighbors and inverse of squared distance
        """

        x_out, y_out, z_out = lat_lon.lon_lat_to_cartesian(lons2d.flatten(), lats2d.flatten())
        dst, ind = self.kdtree.query(zip(x_out, y_out, z_out), k=nneighbours)

        data_in_flat = data_in.flatten()

        inverse_square = 1.0 / dst ** 2
        if len(dst.shape) > 1:
            norm = np.sum(inverse_square, axis=1)
            norm = np.array([norm] * dst.shape[1]).transpose()
            coefs = inverse_square / norm

            data_out_flat = np.sum(coefs * data_in_flat[ind], axis=1)
        elif len(dst.shape) == 1:
            data_out_flat = data_in_flat[ind]
        else:
            raise Exception("Could not find neighbor points")
        return np.reshape(data_out_flat, lons2d.shape)




def create_monthly_means():
    #tmp
    #dm = CRUDataManager()
    #dm.create_monthly_means_file(1901, 2009)

    #pre
    dm = CRUDataManager(path="data/cru_data/CRUTS3.1/cru_ts_3_10.1901.2009.pre.dat.nc", var_name="pre")
    dm.create_monthly_means_file(1901, 2009)


def plot_thawing_index():
    dm = CRUDataManager()
    clim = dm.get_daily_climatology(1981, 2010)
    thi = dm.get_thawing_index_from_climatology(clim)

    plt.pcolormesh(thi.transpose())
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    plot_thawing_index()
    #create_monthly_means()
    #main()
    print "Hello world"
  
