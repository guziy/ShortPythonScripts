import calendar
from collections import OrderedDict
from datetime import datetime
from netCDF4 import Dataset
from matplotlib import gridspec
from matplotlib import cm
from mpl_toolkits.basemap import maskoceans
from matplotlib.colors import BoundaryNorm
from rpn import level_kinds
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

default_season_name_to_months = OrderedDict([
    ("DJF", (12, 1, 2)),
    ("MAM", (3, 4, 5)),
    ("JJA", (6, 7, 8)),
    ("SON", (9, 10, 11))
])

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
                    basemap = rll.get_basemap_object_for_lons_lats(lons2d=lons,
                                                                   lats2d=lats,
                                                                   round=True,
                                                                   resolution="l")
                    r.close()

                    return basemap, lons, lats


def calculate(var_names=None, start_year=None, end_year=None, file_prefix="pm", exp_name=None,
              season_name_to_months=None):
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
        result[vname] = OrderedDict()
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

            ndays_in_season = 0

            for mfp in month_folder_paths:
                for fname in os.listdir(mfp):
                    if not fname.startswith(file_prefix):
                        continue

                    if not fname.endswith("moyenne"):
                        continue

                    the_year, the_month = get_year_and_month_from_name(os.path.basename(mfp))
                    _, days_in_month = calendar.monthrange(the_year, the_month)
                    ndays_in_season += days_in_month
                    fpath = os.path.join(mfp, fname)

                    r = RPN(fpath)
                    for vname in var_names:
                        if vname == "SD":
                            data = r.get_first_record_for_name_and_level(vname,
                                                                         level=1,
                                                                         level_kind=level_kinds.ARBITRARY)
                        elif vname == "AL":
                            data = r.get_first_record_for_name_and_level(vname,
                                                                         level=5,
                                                                         level_kind=level_kinds.ARBITRARY)
                        else:
                            data = r.get_first_record_for_name(vname)

                        vname_to_all_data[vname].append(data * days_in_month)

                    r.close()

            #calculate mean fields for the given year and season
            for vname in var_names:
                #print "{} has {} days".format(season, ndays_in_season)
                result[vname][season].append(np.asarray(vname_to_all_data[vname]).sum(axis=0) / float(ndays_in_season))

        for vname in var_names:
            assert len(result[vname][season]) > 0, "{0}, {1}, {2}, {3}-{4}".format(
                vname, exp_name, season, start_year, end_year)
            assert len(result[vname][season]) == (end_year - start_year + 1)

    return result


def get_clev_by_name(varname, delta=False):
    """

    :rtype : list
    """

    if delta:
        if varname == "I6":
            x = np.arange(-0.1, 0, 0.02)
            return x.tolist() + (-x[::-1]).tolist()
        elif varname == "DN":
            x = np.arange(-100, -25, 25).tolist() + [-15, -10, -5]
            return x + [-xi for xi in x][::-1]
        elif varname == "SD":
            negs = [-50, -30, -25, -20, -15, -10, -5]
            return negs + [-x for x in negs[::-1]]
        elif varname == "AL":
            negs = np.arange(-0.2, 0, 0.02).tolist()
            return negs + [-xi for xi in negs[::-1]]
        elif varname == "I5":
            return np.arange(0, 1000, 100)
    else:

        if varname == "I6":
            return np.arange(0.4, 1, 0.05)
        elif varname == "DN":
            return [0, 10, 20, 50, 100, 200, 400, 600, 700]
        elif varname == "SD":
            return [0, 5, 10, 15, 20, 40, 60, 80, 120, 150]
        elif varname == "AL":
            return np.arange(0., 1, 0.1)
        elif varname == "I5":
            return np.arange(-1000, 1100, 100)


def get_var_label(varname):
    if varname == "I6":
        return "Snow albedo"
    elif varname == "DN":
        return "Snow density"
    elif varname == "SD":
        return "Snow depth"
    elif varname == "AL":
        return "Surface albedo"
    elif varname == "I5":
        return "SWE"


def plot(period_to_data_future=None, period_to_data_current=None, exp_name=None,
         basemap=None, lons=None, lats=None):
    img_folder = "../images/{0}".format(exp_name)
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)

    x, y = basemap(lons, lats)
    font_size = 25
    fig_dpi = 300
    fig_size = (25, 7)
    snowdepth_limit = 1  # everything that is less than 1 cm is ignored

    lons[lons > 180] -= 360
    the_mask = maskoceans(lons, lats, np.zeros(x.shape))

    current_data = None

    period_to_season_to_snow_mass = {}
    period_to_season_to_snow_depth = {}

    #plot seasonal means for the current climate
    for p, data in period_to_data_current.iteritems():
        current_data = data

        #save snow mass and snow depth in order to later calculate snow density
        period_to_season_to_snow_mass[p] = {}
        period_to_season_to_snow_depth[p] = {}

        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=fig_size)
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            clevs = get_clev_by_name(vname, delta=False)
            norm = BoundaryNorm(boundaries=clevs, ncolors=len(clevs) - 1)
            cmap = cm.get_cmap("jet", len(clevs) - 1)

            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0)
                to_plot = np.ma.masked_where(the_mask.mask, to_plot)


                snow_depth = np.asarray(data["SD"][season]).mean(axis=0)
                to_plot = np.ma.masked_where(snow_depth < snowdepth_limit, to_plot)

                if vname == "I5":
                    period_to_season_to_snow_depth[p][season] = to_plot

                if vname == "SD":
                    period_to_season_to_snow_depth[p][season] = to_plot



                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season + "\n", fontsize=font_size)
                im = basemap.pcolormesh(x, y, to_plot, vmin=clevs[0], vmax=clevs[-1],
                                        norm=norm, ax=ax, cmap=cmap)

                #im = basemap.contourf(x, y, to_plot, levels=clevs, norm=norm, ax=ax, extend="max")

                if col == 0:
                    ax.set_ylabel(get_var_label(varname=vname), fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                basemap.drawmapboundary(fill_color="0.75")
                col += 1

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax, extend="max")
            cb.ax.tick_params(labelsize=font_size)
            imname = "{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)

    #future climate
    for p, data in period_to_data_future.iteritems():
        print p, exp_name
        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=fig_size)
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            clevs = get_clev_by_name(vname, delta=False)
            norm = BoundaryNorm(boundaries=clevs, ncolors=len(clevs) - 1)
            cmap = cm.get_cmap("jet", len(clevs) - 1)

            #save snow mass and snow depth in order to later calculate snow density
            period_to_season_to_snow_mass[p] = {}
            period_to_season_to_snow_depth[p] = {}


            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0)
                to_plot = np.ma.masked_where(the_mask.mask, to_plot)

                snow_depth = np.asarray(data["SD"][season]).mean(axis=0)
                to_plot = np.ma.masked_where(snow_depth < snowdepth_limit, to_plot)

                if vname == "I5":
                    period_to_season_to_snow_depth[p][season] = to_plot

                if vname == "SD":
                    period_to_season_to_snow_depth[p][season] = to_plot

                print to_plot.shape, season
                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season + "\n", fontsize=font_size)
                #im = basemap.contourf(x, y, to_plot, levels=clevs, ax=ax, extend="max", norm=norm)
                im = basemap.pcolormesh(x, y, to_plot, vmin=clevs[0], vmax=clevs[-1],
                                        norm=norm, ax=ax, cmap=cmap)

                if col == 0:
                    ax.set_ylabel(get_var_label(varname=vname), fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                basemap.drawmapboundary(fill_color="0.75")
                col += 1

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax, extend="max")
            cb.ax.tick_params(labelsize=font_size)
            imname = "{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)

    #climate change plots

    for p, data in period_to_data_future.iteritems():
        for vname, season_to_data in data.iteritems():

            fig = plt.figure(figsize=fig_size)
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            clevs = get_clev_by_name(vname, delta=True)
            norm = BoundaryNorm(boundaries=clevs, ncolors=len(clevs) - 1)
            diff_cmap = cm.get_cmap("bwr", len(clevs) - 1)
            for season, yearly_fields in season_to_data.iteritems():
                to_plot = np.asarray(yearly_fields).mean(axis=0) - np.asarray(current_data[vname][season]).mean(axis=0)
                to_plot = np.ma.masked_where(the_mask.mask, to_plot)

                snow_depth = np.asarray(data["SD"][season]).mean(axis=0)
                to_plot = np.ma.masked_where(snow_depth < snowdepth_limit, to_plot)

                ax = fig.add_subplot(gs[0, col])
                ax.set_title(season + "\n", fontsize=font_size)
                #im = basemap.contourf(x, y, to_plot, levels=clevs, ax=ax, cmap=diff_cmap,
                #                      extend="both", norm=norm)

                im = basemap.pcolormesh(x, y, to_plot, vmin=clevs[0], vmax=clevs[-1],
                                        norm=norm, ax=ax, cmap=diff_cmap)

                if col == 0:
                    ax.set_ylabel(get_var_label(varname=vname), fontdict={"size": font_size})
                basemap.drawcoastlines(ax=ax)
                col += 1
                basemap.drawmapboundary(fill_color="0.75")

            cax = fig.add_subplot(gs[0, col])

            cb = plt.colorbar(im, cax=cax, ticks=clevs, extend="both")
            cb.ax.tick_params(labelsize=font_size)
            imname = "CC_{0}_{1}_{2}_{3}-{4}.png".format(exp_name, vname, "-".join(season_to_data.keys()), p[0], p[1])
            fig.savefig(os.path.join(img_folder, imname), dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)

        #calculate snow density
        if "SD" in current_data.keys() and "I5" in current_data.keys():
            pass


def calculate_and_plot():
    #var_names = ["DN", "I6", "SD"]
    var_names = ["SD", "AL", "I5"]
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

