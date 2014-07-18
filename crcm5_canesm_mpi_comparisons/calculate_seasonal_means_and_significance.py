import calendar
from collections import OrderedDict
from datetime import datetime
from netCDF4 import Dataset
from matplotlib import gridspec
from matplotlib import cm
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
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

fig_size = (25, 7)


# convenience
period_to_label = {
    (1981, 2010): "C",
    (2041, 2070): "F1",
    (2071, 2100): "F2"
}

CANESM_DRIVEN = "CanESM driven"
CANESM_DRIVEN_DATA_PATH = "/b10_fs1/winger/Arctic/OMSC26_Can_long_new_v01/Diagnostics_free/"

MPI_DRIVEN = "MPI driven"
MPI_DRIVEN_DATA_PATH = "/b10_fs1/winger/Arctic/OMSC26_MPI_long_new_v01/Diagnostics/"

ERA40_DRIVEN = "ERA40 driven"
ERA40_DRIVEN_DATA_PATH = "/b10_fs1/winger/Arctic/Arctic_0.5deg_OMSC_26L_ERA40I"

EXP_to_PATH = {
    CANESM_DRIVEN: CANESM_DRIVEN_DATA_PATH,
    MPI_DRIVEN: MPI_DRIVEN_DATA_PATH
}

default_season_name_to_months = OrderedDict([
    ("DJF", (12, 1, 2)),
    ("MAM", (3, 4, 5)),
    ("JJA", (6, 7, 8)),
    ("SON", (9, 10, 11))
])

nc_cache_folder = "../nc_cache"


def get_land_sea_glaciers_mask_from_geophysics_file(
        path="/b10_fs1/winger/Arctic/OMSC26_Can_long_new_v01/Geophys/land_sea_glacier_mask_free"):
    r = RPN(path)
    mask = r.get_first_record_for_name("FMSK") < 0.5
    r.close()
    return mask


def get_basemap_and_coordinates_from_any_file_in(folder=CANESM_DRIVEN_DATA_PATH, prefix="pm"):
    for monthname in os.listdir(folder):
        monthpath = os.path.join(folder, monthname)
        if os.path.isdir(monthpath):
            for fname in os.listdir(monthpath):
                if fname.startswith(prefix) and fname.lower().endswith("moyenne"):
                    fpath = os.path.join(monthpath, fname)
                    r = RPN(fpath)

                    vnames = r.get_list_of_varnames()
                    any_name = [v for v in vnames if v not in [">>", "^^", "HY"]][0]
                    r.get_first_record_for_name(any_name)

                    #rll = RotatedLatLon(**r.get_proj_parameters_for_the_last_read_rec())
                    lons, lats = r.get_longitudes_and_latitudes_for_the_last_read_rec()
                    #basemap = rll.get_basemap_object_for_lons_lats(lons2d=lons,
                    #                                               lats2d=lats,
                    #                                               round=True,
                    #                                               resolution="l")

                    basemap = Basemap(projection="npstere", lon_0=10, boundinglat=45, round=True)

                    r.close()

                    return basemap, lons, lats


def calculate(var_names=None, start_year=None, end_year=None, file_prefix="pm",
              exp_name=None, exp_data_path=None,
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

    data_folder = EXP_to_PATH[exp_name] if exp_data_path is None else exp_data_path

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

                    if file_prefix is not None and not fname.startswith(file_prefix):
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
                        elif vname == "SN":
                            data = r.get_all_time_records_for_name(vname)
                        elif vname == "TT":
                            r1 = RPN(os.path.join(mfp, "dm" + fname[2:]))
                            data = r1.get_first_record_for_name_and_level(vname, level=1,
                                                                          level_kind=level_kinds.HYBRID)
                            r1.close()

                        else:
                            data = r.get_first_record_for_name(vname)

                        # several fields with the same name and level in the same file (but of different shapes)
                        if isinstance(data, dict):
                            selected_field = None
                            for d, field in data.iteritems():
                                if field.shape == (172, 160):
                                    selected_field = field
                                    break

                            data = selected_field

                        if vname in ["SN", "PR"]:
                            data *= 1000 * 60 * 60 * 24

                        if data.shape == (198, 186):
                            data = data[13:-13, 13:-13]

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
            negs = np.arange(-0.2, 0, 0.04).tolist()
            return negs + [-xi for xi in negs[::-1]]
        elif varname == "I5":
            x = np.arange(-200, 0, 50).tolist() + [-25, -5]
            return x + [-xi for xi in reversed(x)]
        elif varname == "SN":
            x = np.arange(-1, 0, 0.25).tolist() + [-0.05, ]
            return x + [-xi for xi in reversed(x)]
        elif varname == "TT":
            x = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
            return x
        elif varname == "PR":
            x = [-3, -2.5, -1.5, -0.8, -0.4, -0.2]
            return x + [-xi for xi in reversed(x)]
        elif varname == "AH":
            x = [-100, -80, -50, -10, -5, -1]
            return x + [-xi for xi in reversed(x)]
        elif varname == "AV":
            x = [-100, -80, -50, -10, -5, -1]
            return x + [-xi for xi in reversed(x)]


    else:
        if varname == "I6":
            return np.arange(0.4, 1, 0.05)
        elif varname == "DN":
            return np.arange(50, 700, 50)
        elif varname == "SD":
            return [0, 5, 10, 15, 20, 40, 60, 80, 120, 150]
        elif varname == "AL":
            return np.arange(0., 1, 0.1)
        elif varname == "I5":
            return np.arange(0, 300, 25)
        elif varname == "SN":
            return np.arange(0, 2.8, 0.2)
        elif varname == "TT":
            return [-25, -15, -10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10, 15, 25]
        elif varname == "PR":
            return [0, 0.2, 0.5, 1, 2, 5, 8, 10]
        elif varname in ["AH"]:
            return [-20, -10, 0, 10, 30, 50, 100, 150]
        elif varname == "AV":
            return [0, 10, 30, 50, 80, 100, 150]


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
    elif varname == "SN":
        return "Snow fall"
    elif varname == "TT":
        return "2m Air temperature"
    elif varname == "AH":
        return "Sensible heat flux"
    elif varname == "AV":
        return "Latent heat flux"
    else:
        return varname


def plot(period_to_data_future=None, period_to_data_current=None, exp_name=None,
         basemap=None, lons=None, lats=None):
    img_folder = "../images/{0}".format(exp_name)
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)

    x, y = basemap(lons, lats)

    fig_dpi = 100
    snowdepth_limit = 1  # everything that is less than 1 cm is ignored

    #lons[lons > 180] -= 360
    #the_mask = maskoceans(lons, lats, np.zeros(x.shape))

    current_data = None
    current_period = None

    period_to_season_to_snow_mass = {}
    period_to_season_to_snow_depth = {}

    land_sea_glaciers_mask = get_land_sea_glaciers_mask_from_geophysics_file()
    #	path="/b10_fs1/winger/Arctic/OMSC26_Can_long_new_v01/Geophys/land_sea_mask_free")

    #plot seasonal means for the current climate
    for p, data in period_to_data_current.iteritems():
        current_data = data
        current_period = p

        #save snow mass and snow depth in order to later calculate snow density
        period_to_season_to_snow_mass[p] = OrderedDict()
        period_to_season_to_snow_depth[p] = OrderedDict()

        print data.keys()

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
                to_plot = np.ma.masked_where(land_sea_glaciers_mask, to_plot)

                snow_depth = np.asarray(data["SD"][season]).mean(axis=0)
                if vname not in ["TT", "AH", "AV"]:
                    to_plot = np.ma.masked_where(snow_depth < snowdepth_limit, to_plot)

                if vname == "I5":
                    period_to_season_to_snow_mass[p][season] = to_plot

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
        print data.keys()
        #save snow mass and snow depth in order to later calculate snow density
        period_to_season_to_snow_mass[p] = OrderedDict()
        period_to_season_to_snow_depth[p] = OrderedDict()

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
                to_plot = np.ma.masked_where(land_sea_glaciers_mask, to_plot)

                if vname not in ["TT", "AH", "AV"]:
                    snow_depth = np.asarray(data["SD"][season]).mean(axis=0)
                    to_plot = np.ma.masked_where(snow_depth < snowdepth_limit, to_plot)

                if vname == "I5":
                    period_to_season_to_snow_mass[p][season] = to_plot

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
                to_plot = np.ma.masked_where(land_sea_glaciers_mask, to_plot)

                if vname not in ["TT", "AH", "AV"]:
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

    for p in period_to_season_to_snow_mass:
        print "Snow mass", p, period_to_season_to_snow_mass[p].keys()
        print "Snow depth", p, period_to_season_to_snow_depth[p].keys()


    #calculate snow density
    if "SD" in current_data.keys() and "I5" in current_data.keys():
        ref_snow_depth = period_to_season_to_snow_depth[current_period]
        ref_snow_mass = period_to_season_to_snow_mass[current_period]
        vname = "DN"
        #plot values
        for p in period_to_season_to_snow_depth:
            fig = plt.figure(figsize=fig_size)
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            clevs = get_clev_by_name(vname, delta=False)
            norm = BoundaryNorm(boundaries=clevs, ncolors=len(clevs) - 1)
            cmap = cm.get_cmap("jet", len(clevs) - 1)

            season_to_snowdepth = period_to_season_to_snow_depth[p]
            season_to_snowmass = period_to_season_to_snow_mass[p]

            for season, snow_depth in season_to_snowdepth.iteritems():

                print season_to_snowmass.keys(), season_to_snowdepth.keys()
                #calculate snow density
                to_plot = season_to_snowmass[season] / (1e-2 * snow_depth)

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
            imname = "{0}_{1}_{2}_{3}-{4}.png".format(exp_name,
                                                      vname + "_calculated", "-".join(season_to_data.keys()),
                                                      p[0], p[1])

            fig.savefig(os.path.join(img_folder, imname), dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)

        #plot climate change
        for p in period_to_season_to_snow_depth:

            season_to_snowdepth = period_to_season_to_snow_depth[p]
            season_to_snowmass = period_to_season_to_snow_mass[p]

            #  No need to calculate climate change in cur climate with respect to the current climate
            if season_to_snowdepth is ref_snow_depth:
                continue

            fig = plt.figure(figsize=fig_size)
            gs = gridspec.GridSpec(nrows=1, ncols=len(season_to_data) + 1,
                                   width_ratios=len(season_to_data) * [1.0, ] + [0.05, ])  # +1 for colorbar
            col = 0
            im = None
            clevs = get_clev_by_name(vname, delta=True)
            norm = BoundaryNorm(boundaries=clevs, ncolors=len(clevs) - 1)
            cmap = cm.get_cmap("bwr", len(clevs) - 1)

            for season, snow_depth in season_to_snowdepth.iteritems():
                #calculate snow density
                to_plot = season_to_snowmass[season] / (1e-2 * snow_depth) - ref_snow_mass[season] / (
                    1e-2 * ref_snow_depth[season])

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

            cb = plt.colorbar(im, cax=cax, extend="both", ticks=clevs)
            cb.ax.tick_params(labelsize=font_size)
            imname = "CC_{0}_{1}_{2}_{3}-{4}.png".format(exp_name,
                                                         vname + "_calculated", "-".join(season_to_data.keys()),
                                                         p[0], p[1])

            fig.savefig(os.path.join(img_folder, imname), dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)


def calculate_snowdepth_max_for_seasons(season_to_months=None, exp_name=None, start_year=None, end_year=None,
                                        file_name_prefix="pm", file_name_suffix="moyenne"):
    """
    returns a dictionary {year: {season: max snow depth during the season}}
    :param season_to_months:
    :param exp_name:
    :param start_year:
    :param end_year:
    :param file_name_prefix:
    :param file_name_suffix:
    :return:
    """
    if season_to_months is None:
        season_to_months = OrderedDict(
            ("NDJFM", (11, 12, 1, 2, 3),)
        )

    data_path = EXP_to_PATH[exp_name]
    monthly_folder_names = [mfn for mfn in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, mfn))]
    import pandas as pd

    mf_names = pd.Series(index=monthly_folder_names, data=monthly_folder_names)

    def grouper(mfname):
        the_year, the_month = get_year_and_month_from_name(mfname)
        the_season = None

        if start_year <= the_year <= end_year:
            for a_season, months in season_to_months.iteritems():
                if the_month in months:
                    the_season = a_season
                    break

        return the_year, the_season

    mf_names_groups = mf_names.groupby(by=grouper)

    year_to_season_to_max_snow_depth = OrderedDict()
    for year in range(start_year, end_year + 1):
        year_to_season_to_max_snow_depth[year] = OrderedDict()

    #Calculate maximum snow depth field for each season of each month
    for key, val in mf_names_groups:
        year, season = key

        if season is None:
            continue

        fields_for_year = []

        current_monthly_folders = [os.path.join(data_path, mfn) for mfn in val]
        for mfp in current_monthly_folders:
            for fname in os.listdir(mfp):

                if not fname.startswith(file_name_prefix):
                    continue

                if not fname.endswith(file_name_suffix):
                    continue

                fpath = os.path.join(mfp, fname)

                r = RPN(fpath)
                fields_for_year.append(r.get_first_record_for_name_and_level(varname="SD", level=1))
                r.close()

        year_to_season_to_max_snow_depth[year][season] = np.max(fields_for_year, axis=0)

    return year_to_season_to_max_snow_depth


def plot_max_snow_panel(year_to_season_to_max_snow_depth=None, exp_name=None,
                        basemap=None, lons=None, lats=None, period_current=None, periods_future=None):
    img_folder = "../images/{0}".format(exp_name)
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)

    x, y = basemap(lons, lats)

    fig_dpi = 100
    fig_size_sdmax = (25, 14)
    snowdepth_limit = 1  # everything that is less than 1 cm is ignored

    current_data = None

    land_sea_glaciers_mask = get_land_sea_glaciers_mask_from_geophysics_file()

    seasons = year_to_season_to_max_snow_depth.items()[0][1].keys()
    n_seasons = len(seasons)

    clevs_diff = get_clev_by_name("SD", delta=True)
    norm_diff = BoundaryNorm(boundaries=clevs_diff, ncolors=len(clevs_diff) - 1)
    cmap_diff = cm.get_cmap("bwr", len(clevs_diff) - 1)

    clevs_val = get_clev_by_name("SD", delta=False)
    norm_val = BoundaryNorm(boundaries=clevs_val, ncolors=len(clevs_val) - 1)
    cmap_val = cm.get_cmap("jet", len(clevs_val) - 1)

    fig = plt.figure(figsize=fig_size_sdmax)
    gs = gridspec.GridSpec(nrows=n_seasons + 1, ncols=len(periods_future) + 2, width_ratios=[1, 1, 1, 0.05])

    #plot current climate
    current_season_to_mean_sdmax = OrderedDict()
    sel_data = [year_to_season_to_max_snow_depth[year] for year in range(period_current[0], period_current[1] + 1)]

    row = 0
    for the_season in seasons:
        to_plot = np.asarray([season_to_field[the_season] for season_to_field in sel_data]).mean(axis=0)

        to_plot = np.ma.masked_where(land_sea_glaciers_mask | (to_plot < snowdepth_limit), to_plot)

        current_season_to_mean_sdmax[the_season] = to_plot

        ax = fig.add_subplot(gs[row, 0])
        if row == 0:
            ax.set_title(period_to_label[period_current])
        basemap.pcolormesh(x, y, current_season_to_mean_sdmax[the_season], vmin=clevs_val[0], vmax=clevs_val[-1],
                           norm=norm_val, cmap=cmap_val, ax=ax)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")
        row += 2

    #plot future climate
    period_to_season_to_sdmax = OrderedDict()
    col = 1
    im = None
    for p in periods_future:
        future_season_to_mean_sdmax = OrderedDict()
        period_to_season_to_sdmax[p] = future_season_to_mean_sdmax
        sel_data = [year_to_season_to_max_snow_depth[year] for year in range(p[0], p[1] + 1)]

        row = 0
        for the_season in seasons:
            to_plot = np.asarray([season_to_field[the_season] for season_to_field in sel_data]).mean(axis=0)

            ax = fig.add_subplot(gs[row, col])
            if row == 0:
                ax.set_title(period_to_label[p])
            to_plot = np.ma.masked_where(land_sea_glaciers_mask | (to_plot < 1), to_plot)
            future_season_to_mean_sdmax[the_season] = to_plot
            im = basemap.pcolormesh(x, y, to_plot,
                                    vmin=clevs_val[0], vmax=clevs_val[-1],
                                    norm=norm_val, cmap=cmap_val, ax=ax)
            basemap.drawcoastlines()
            basemap.drawmapboundary(fill_color="0.75")
            row += 2
        col += 1

    cax = fig.add_subplot(gs[0, col])
    plt.colorbar(im, cax=cax, extend="max")

    #plot climate change
    col = 1
    for p in periods_future:
        row = 1
        for season in seasons:

            cc = period_to_season_to_sdmax[p][season] - current_season_to_mean_sdmax[season]
            print row, col
            ax = fig.add_subplot(gs[row, col])
            im = basemap.pcolormesh(x, y, cc, vmin=clevs_diff[0], vmax=clevs_diff[-1],
                                    norm=norm_diff, cmap=cmap_diff, ax=ax)
            basemap.drawcoastlines()
            basemap.drawmapboundary(fill_color="0.75")
            if row == 1:
                ax.set_title("{}-{}".format(period_to_label[p], period_to_label[period_current]))

            row += 2
        col += 1

    cax = fig.add_subplot(gs[1, col])
    plt.colorbar(im, cax=cax, extend="both")

    figname_elts = ("-".join(seasons), exp_name) + period_current + tuple(np.ravel(periods_future).tolist()) + ("SD", )
    fig.savefig(os.path.join(img_folder, "{}_{}_({}-{})_({}-{})_({}-{})_{}max.png".format(*figname_elts)),
                dpi=fig_dpi)


def check_calculate_and_plot_snowdepth_max_for_season():
    start_year = 1985
    end_year = 2000
    calculate_snowdepth_max_for_seasons(exp_name=CANESM_DRIVEN, start_year=start_year, end_year=end_year)


def calculate_and_plot():
    #var_names = ["DN", "I6", "SD"]
    #var_names = ["SD", "AL", "I5", "SN", "DN", "I6"]
    #var_names = ["SD", "SN"]
    file_prefix = "pm"

    var_names = ["SD", "TT", "AH", "AV"]

    period_c = (1981, 2010)
    period_f1 = (2041, 2070)
    period_f2 = (2071, 2100)

    season_to_months_for_max_snow_depth = \
        OrderedDict(
            [("NDJFM", [11, 12, 1, 2, 3])])

    season_to_months = default_season_name_to_months

    periods_future = [period_f1, period_f2]

    #coordinates
    basemap, lons, lats = get_basemap_and_coordinates_from_any_file_in(folder=CANESM_DRIVEN_DATA_PATH)

    experiments = [MPI_DRIVEN, CANESM_DRIVEN]
    for exp in experiments:
        period_to_data_future = {}
        data_current = calculate(var_names, start_year=period_c[0], end_year=period_c[1],
                                 file_prefix=file_prefix, exp_name=exp,
                                 season_name_to_months=season_to_months)

        for p in periods_future:
            period_to_data_future[p] = calculate(var_names, start_year=p[0], end_year=p[1],
                                                 file_prefix=file_prefix, exp_name=exp,
                                                 season_name_to_months=season_to_months)

        period_to_data_current = {
            period_c: data_current
        }
        plot(period_to_data_future=period_to_data_future, period_to_data_current=period_to_data_current,
             exp_name=exp, basemap=basemap, lons=lons, lats=lats)

    #Plot max snow depth
    for exp in experiments:
        year_to_season_to_data = calculate_snowdepth_max_for_seasons(
            season_to_months=season_to_months_for_max_snow_depth,
            exp_name=exp,
            start_year=period_c[0], end_year=periods_future[-1][1])

        plot_max_snow_panel(year_to_season_to_max_snow_depth=year_to_season_to_data, exp_name=exp, basemap=basemap,
                            lons=lons, lats=lats, period_current=period_c, periods_future=periods_future)


if __name__ == "__main__":
    assert os.path.isdir(MPI_DRIVEN_DATA_PATH) and os.path.isdir(MPI_DRIVEN_DATA_PATH)
    import matplotlib.pyplot as plt

    font_size = 25
    plt.rcParams.update({
        'axes.labelsize': font_size,
        'font.size': font_size,
        'text.fontsize': font_size,
        'legend.fontsize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'figure.figsize': fig_size,
        "axes.titlesize": font_size
    })
    calculate_and_plot()
    #check_calculate_and_plot_snowdepth_max_for_season()

