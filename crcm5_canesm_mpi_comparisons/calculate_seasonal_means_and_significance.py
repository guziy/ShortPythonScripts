from collections import OrderedDict
from datetime import datetime
from netCDF4 import Dataset
from matplotlib import gridspec
from matplotlib import cm
from rpn.rpn import RPN
import numpy as np
from crcm5_canesm_mpi_comparisons.util import get_year_and_month_from_name
from domains.rotated_lat_lon import RotatedLatLon
import matplotlib.pyplot as plt


__author__ = 'huziy'

"""
Calculate seasonal means and save them to the cached files for plotting later, plotting procedures
are also in this module

Uses monthly means as an input
"""

import os

CANESM_DRIVEN = "CanESM driven"
CANESM_DATA_PATH = "/b10_fs1/winger/Arctic/OMSC26_Can_long_new_v01/Diagnostics_free/"

MPI_DRIVEN = "MPI driven"
MPI_DATA_PATH = "/b10_fs1/winger/Arctic/OMSC26_MPI_long_new_v01/Diagnostics/"

EXP_to_PATH = {
    CANESM_DRIVEN: CANESM_DATA_PATH,
    MPI_DRIVEN: MPI_DATA_PATH
}

default_season_name_to_months = OrderedDict({
    "DJF": (12, 1, 2),
    "MAM": (3, 4, 5),
    "JJA": (6, 7, 8),
    "SON": (9, 10, 11)
})

nc_cache_folder = "../nc_cache"


def get_basemap_and_coordinates_from_any_file_in(folder=CANESM_DATA_PATH, prefix="pm"):
    for monthname in os.listdir(folder):
        monthpath = os.path.join(CANESM_DATA_PATH, monthname)
        if os.path.isdir(monthpath):
            for fname in os.listdir(monthpath):
                if fname.startswith(prefix) and fname.lower().endswith("moyenne"):
                    fpath = os.path.join(monthpath, fname)
                    r = RPN(fpath)

                    vnames = r.get_list_of_varnames()
                    any_name = [v for v in vnames if v not in [">>", "^^", "HY"]][0]
                    r.get_first_record_for_name(any_name)
                    rll = RotatedLatLon(**r.get_proj_parameters_for_the_last_read_rec())
                    lons, lats = r.get_longitudes_and_latitudes_for_the_last_read_rec()
                    basemap = rll.get_basemap_object_for_lons_lats(lons2d=lons, lats2d=lats)
                    r.close()

                    return basemap, lons, lats


def calculate(var_names=None, start_year=None, end_year=None, file_prefix="pm", exp_name=None, season_name_to_months=None):
    """
    Calculate seasonal means and store it to the cache files

    return {varname: {season_name: 3d field with dimensions (year, x, y)}}


    :param var_names:
    :param start_year:
    :param end_year:
    """

    if not os.path.isdir(nc_cache_folder):
        os.mkdir(nc_cache_folder)

    if season_name_to_months is None:
        season_name_to_months = default_season_name_to_months

    data_folder = EXP_to_PATH[exp_name]

    result = {}

    for vname in var_names:
        result[vname] = OrderedDict({})
        for season, month_list in season_name_to_months.iteritems():
            result[vname][season] = []

    for season, month_list in season_name_to_months.iteritems():
        for year in range(start_year, end_year + 1):

            month_folder_names = [mfn for mfn in os.listdir(data_folder)
                                  if get_year_and_month_from_name(fname=mfn)[0] == year and
                                     get_year_and_month_from_name(fname=mfn)[1] in month_list]

            month_folder_paths = [os.path.join(data_folder, mfn) for mfn in month_folder_names]

            vname_to_all_data = {}
            for vname in var_names:
                vname_to_all_data[vname] = []

            for mfp in month_folder_paths:
                for fname in os.listdir(mfp):
                    if not fname.startswith(file_prefix):
                        continue
                    fpath = os.path.join(mfp, fname)

                    r = RPN(fpath)
                    for vname in var_names:
                        vname_to_all_data[vname].append(r.get_first_record_for_name(vname))
                    r.close()

            #calculate mean fields for the given year and season
            for vname in var_names:
                result[vname][season].append(np.asarray(vname_to_all_data[vname]).mean(axis=0))

        for vname in var_names:
            assert len(result[vname][season]) > 0, "{0}, {1}, {2}, {3}-{4}".format(
                vname, exp_name, season, start_year, end_year)

    return result


def get_clev_by_name(varname, delta=False):
    if varname == "I6" and not delta:
        return np.arange(0, 1, 0.05)
    elif varname == "DN" and not delta:
        return np.arange(0, 600, 50)
    elif varname == "I6" and delta:
        return np.arange(-1, 1.1, 0.1)
    elif varname == "DN":
        return np.arange(-600, 650, 50)


def plot(period_to_data_future=None, period_to_data_current=None, exp_name=None,
         basemap=None, lons=None, lats=None):
    img_folder = "../images/{0}".format(exp_name)
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)

    x, y = basemap(lons, lats)
    font_size = 25

    current_data = None
    #plot seasonal means for the current climate
    for p, data in period_to_data_current.iteritems():
        current_data = data
        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=(20, 5))
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0)
                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season, fontsize=font_size)
                im = basemap.contourf(x, y, to_plot, levels=get_clev_by_name(vname), ax=ax)
                ax.set_ylabel(vname, fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                col += 1

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=font_size)
            imname = "{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=800)

    #future climate
    for p, data in period_to_data_future.iteritems():
        print p, exp_name
        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=(20, 5))
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0)
                print to_plot.shape, season
                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season, fontsize=font_size)
                im = basemap.contourf(x, y, to_plot, levels=get_clev_by_name(vname), ax=ax)
                ax.set_ylabel(vname, fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                col += 1

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=font_size)
            imname = "{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=800)

    #climate change plots
    diff_cmap = cm.get_cmap("RdBu_r")
    for p, data in period_to_data_future.iteritems():
        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=(20, 5))
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0) - np.asarray(current_data[vname][season]).mean(axis=0)
                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season, fontsize=font_size)
                im = basemap.contourf(x, y, to_plot, levels=get_clev_by_name(vname), ax=ax, cmap=diff_cmap)
                ax.set_ylabel(vname, fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                col += 1

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=font_size)
            imname = "CC_{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=800)


def calculate_and_plot():
    var_names = ["DN", "I6"]
    file_prefix = "pm"

    period_c = (1981, 2010)
    period_f1 = (2041, 2070)
    period_f2 = (2071, 2100)

    periods_future = [period_f1, period_f2]

    #coordinates
    basemap, lons, lats = get_basemap_and_coordinates_from_any_file_in(folder=CANESM_DATA_PATH)

    experiments = [MPI_DRIVEN, CANESM_DRIVEN]
    for exp in experiments:
        period_to_data_future = {}
        data_current = calculate(var_names, start_year=period_c[0], end_year=period_c[1],
                                 file_prefix=file_prefix, exp_name=exp)

        for p in periods_future:
            period_to_data_future[p] = calculate(var_names, start_year=p[0], end_year=p[1],
                                                 file_prefix=file_prefix, exp_name=exp)
        period_to_data_current = {
            period_c: data_current
        }
        plot(period_to_data_future=period_to_data_future, period_to_data_current=period_to_data_current,
             exp_name=exp, basemap=basemap, lons=lons, lats=lats)


if __name__ == "__main__":
    assert os.path.isdir(MPI_DATA_PATH) and os.path.isdir(MPI_DATA_PATH)
    calculate_and_plot()

