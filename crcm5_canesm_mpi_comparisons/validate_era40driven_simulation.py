from collections import OrderedDict
import os
from matplotlib import cm
from matplotlib.colors import BoundaryNorm
from matplotlib.gridspec import GridSpec
from crcm5_canesm_mpi_comparisons.cmc_swe_snow import CmcSweSnowManager
from cru.temperature import CRUDataManager

__author__ = 'huziy'

import calculate_seasonal_means_and_significance as model_data_calculator
import numpy as np
import matplotlib.pyplot as plt

ERA40_DRIVEN = "ERA40 driven"
ERA40_DRIVEN_DATA_PATH = "/b10_fs1/winger/Arctic/Arctic_0.5deg_OMSC_26L_ERA40I/Diagnostics"

EXP_to_PATH = {
    ERA40_DRIVEN: ERA40_DRIVEN_DATA_PATH
}

default_season_name_to_months = OrderedDict([
    ("DJF", (12, 1, 2)),
    ("MAM", (3, 4, 5)),
    ("JJA", (6, 7, 8)),
    ("SON", (9, 10, 11))
])


def validate_albedo(start_year=None, end_year=None):
    varnames = ["AL"]
    model_albedo = model_data_calculator.calculate(var_names=varnames,
                                                   start_year=start_year, end_year=end_year,
                                                   exp_name=ERA40_DRIVEN, exp_data_path=ERA40_DRIVEN_DATA_PATH,
                                                   season_name_to_months=default_season_name_to_months)

    #calculate climatology from model data
    model_albedo_clim = OrderedDict()
    for vname, season_to_data_for_years in model_albedo.iteritems():
        for season in season_to_data_for_years:
            model_albedo_clim[season] = np.mean(season_to_data_for_years[season])


            #get and calculate climatology from modis data
            #TODO:


def validate_swe_and_snowdepth(start_year=None, end_year=None, basemap=None, x=None, y=None, gl_land_sea_mask=None):
    """
    Compares climatological means using CMC analysis data
    """
    lons_model, lats_model = basemap(x, y, inverse=True)

    varnames = ["SD"]
    model_data = model_data_calculator.calculate(var_names=varnames,
                                                 start_year=start_year, end_year=end_year,
                                                 exp_name=ERA40_DRIVEN, exp_data_path=ERA40_DRIVEN_DATA_PATH,
                                                 season_name_to_months=default_season_name_to_months)

    #calculate climatology from model data
    model_clim = OrderedDict()
    for vname, season_to_data_for_years in model_data.iteritems():
        model_clim[vname] = OrderedDict()
        for season in season_to_data_for_years:
            model_clim[vname][season] = np.mean(season_to_data_for_years[season], axis=0)

    path_snow_depth = "/b2_fs2/huziy/PythonProjects/DATA/CMC/nsidc0447_CMC_snow_depth_v01/Monthly_Snow_Depth_Files"
    path_swe = "/b2_fs2/huziy/PythonProjects/DATA/CMC/nsidc0447_CMC_snow_depth_v01/Monthly_SWE_Estimates"

    cmc = CmcSweSnowManager(path_snow_depth=path_snow_depth, path_swe=path_swe)


    #snow depth mod, obs climatology and biases
    fig = plt.figure()
    nseasons = len(default_season_name_to_months)
    gs = GridSpec(3, nseasons + 1, width_ratios=[1.0, ] * nseasons + [0.05, ])

    clevs_snowdepth = model_data_calculator.get_clev_by_name("SD")
    cmap_snowdepth = cm.get_cmap("jet", len(clevs_snowdepth) - 1)
    norm_snowdepth = BoundaryNorm(clevs_snowdepth, len(clevs_snowdepth) - 1)

    clevs_snowdepth_diff = model_data_calculator.get_clev_by_name("SD", delta=True)
    cmap_snowdepth_diff = cm.get_cmap("bwr", len(clevs_snowdepth_diff) - 1)
    norm_snowdepth_diff = BoundaryNorm(clevs_snowdepth_diff, len(clevs_snowdepth_diff) - 1)

    im = None

    #plot model climatology
    col = 0
    for season, data in model_clim["SD"].iteritems():
        ax = fig.add_subplot(gs[0, col])

        if col == 0:
            ax.set_ylabel("Model (ERA40-driven)")

        ax.set_title(season)

        to_plot = np.ma.masked_where(gl_land_sea_mask, data)
        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs_snowdepth[0], vmax=clevs_snowdepth[-1], cmap=cmap_snowdepth,
                                norm=norm_snowdepth)
        basemap.drawmapboundary(fill_color="0.75")
        basemap.drawcoastlines()
        col += 1




    #plot obs climatology
    obs_snow_clim = cmc.get_seasonal_mean_clim(start_year=start_year, end_year=end_year,
                                               season_name_to_months=default_season_name_to_months,
                                               var_name="sd")

    col = 0
    bias_clim_snowdepth = OrderedDict()
    for season, data in obs_snow_clim.iteritems():

        ax = fig.add_subplot(gs[1, col])

        if col == 0:
            ax.set_ylabel("Obs(CMC)")

        to_plot = cmc.interpolate_data_to(target_lons=lons_model, target_lats=lats_model, data_field=data,
                                          nneighbours=4)

        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs_snowdepth[0], vmax=clevs_snowdepth[-1], cmap=cmap_snowdepth,
                                norm=norm_snowdepth)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")

        bias_clim_snowdepth[season] = model_clim["SD"][season] - to_plot
        col += 1

    cax = fig.add_subplot(gs[:2, col])
    plt.colorbar(im, ticks=clevs_snowdepth, cax=cax, extend="max")


    #plot mod-obs climatology
    col = 0
    im = None
    for season, data in bias_clim_snowdepth.iteritems():

        ax = fig.add_subplot(gs[2, col])

        if col == 0:
            ax.set_ylabel("Model - Obs")

        im = basemap.pcolormesh(x, y, data, vmin=clevs_snowdepth_diff[0],
                                vmax=clevs_snowdepth_diff[-1], cmap=cmap_snowdepth_diff,
                                norm=norm_snowdepth_diff)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")
        col += 1

    plt.colorbar(im, ticks=clevs_snowdepth_diff, extend="both", cax=fig.add_subplot(gs[2, col]))

    img_folder = "../images/validation/cmc/snow_depth"
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)
    img_path = os.path.join(img_folder, "{}-{}.png".format(start_year, end_year))
    fig.savefig(img_path, dpi=100)
    plt.close(fig)

    ####Plot SWE data
    fig = plt.figure(figsize=(13, 18))

    season_to_months_swe = OrderedDict((
        ("DJF", (1, 2, 12)),
        ("MAM", (3, 4, 5))))

    varnames = ["I5"]
    model_data = model_data_calculator.calculate(var_names=varnames,
                                                 start_year=start_year, end_year=end_year,
                                                 exp_name=ERA40_DRIVEN, exp_data_path=ERA40_DRIVEN_DATA_PATH,
                                                 season_name_to_months=season_to_months_swe)

    #calculate climatology from model data
    model_clim = OrderedDict()
    for vname, season_to_data_for_years in model_data.iteritems():
        model_clim[vname] = OrderedDict()
        for season in season_to_data_for_years:
            model_clim[vname][season] = np.mean(season_to_data_for_years[season], axis=0)

    nseasons = len(season_to_months_swe)
    gs = GridSpec(3, nseasons + 1, width_ratios=[1.0, ] * nseasons + [0.05, ])
    clevs_swe = model_data_calculator.get_clev_by_name("I5")
    cmap_swe = cm.get_cmap("jet", len(clevs_swe) - 1)
    norm_swe = BoundaryNorm(clevs_swe, len(clevs_swe) - 1)

    clevs_swe_diff = model_data_calculator.get_clev_by_name("I5", delta=True)
    cmap_swe_diff = cm.get_cmap("bwr", len(clevs_swe_diff) - 1)
    norm_swe_diff = BoundaryNorm(clevs_swe_diff, len(clevs_swe_diff) - 1)

    im = None

    #plot model climatology
    col = 0
    for season, data in model_clim["I5"].iteritems():
        print col
        ax = fig.add_subplot(gs[0, col])

        if col == 0:
            ax.set_ylabel("Model (ERA40-driven)")

        ax.set_title(season)

        to_plot = np.ma.masked_where(gl_land_sea_mask, data)
        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs_swe[0], vmax=clevs_swe[-1], cmap=cmap_swe,
                                norm=norm_swe)
        basemap.drawmapboundary(fill_color="0.75")
        basemap.drawcoastlines()
        col += 1

    #plot obs climatology
    obs_snow_clim = cmc.get_seasonal_mean_clim(start_year=start_year, end_year=end_year,
                                               season_name_to_months=season_to_months_swe,
                                               var_name="swe")
    col = 0
    bias_clim_swe = OrderedDict()
    for season, data in obs_snow_clim.iteritems():

        ax = fig.add_subplot(gs[1, col])

        if col == 0:
            ax.set_ylabel("Obs (CMC)")

        to_plot = cmc.interpolate_data_to(target_lons=lons_model, target_lats=lats_model, data_field=data,
                                          nneighbours=4)

        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs_swe[0], vmax=clevs_swe[-1], cmap=cmap_swe,
                                norm=norm_swe)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")

        bias_clim_swe[season] = model_clim["I5"][season] - to_plot
        col += 1

    cax = fig.add_subplot(gs[:2, col])
    plt.colorbar(im, ticks=clevs_swe, cax=cax, extend="max")


    #plot mod-obs climatology
    col = 0
    im = None
    for season, data in bias_clim_swe.iteritems():

        ax = fig.add_subplot(gs[2, col])

        if col == 0:
            ax.set_ylabel("Model - Obs")

        im = basemap.pcolormesh(x, y, data, vmin=clevs_swe_diff[0],
                                vmax=clevs_swe_diff[-1], cmap=cmap_swe_diff,
                                norm=norm_swe_diff)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")
        col += 1

    plt.colorbar(im, ticks=clevs_swe_diff, extend="both", cax=fig.add_subplot(gs[2, col]))

    img_folder = "../images/validation/cmc/swe"
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)
    img_path = os.path.join(img_folder, "{}-{}.png".format(start_year, end_year))
    fig.savefig(img_path, dpi=100)
    plt.close(fig)


def validate_with_cru(start_year=None, end_year=None, lons=None, lats=None, basemap=None,
                      gl_land_sea_mask=None, season_to_months=None, model_var_name="TT",
                      cru_var_name="tmp", cru_file=""):
    """
    Using CRU data

    """

    x, y = basemap(lons, lats)

    cru = CRUDataManager(var_name=cru_var_name, path=cru_file)
    seasonal_obs_data = cru.get_seasonal_means(season_name_to_months=default_season_name_to_months,
                                               start_year=start_year,
                                               end_year=end_year)

    varnames = [model_var_name, ]
    model_data = model_data_calculator.calculate(var_names=varnames,
                                                 start_year=start_year, end_year=end_year,
                                                 exp_name=ERA40_DRIVEN, exp_data_path=ERA40_DRIVEN_DATA_PATH,
                                                 season_name_to_months=season_to_months)

    #calculate climatology from model data
    model_clim = OrderedDict()
    for vname, season_to_data_for_years in model_data.iteritems():
        model_clim[vname] = OrderedDict()
        for season in season_to_data_for_years:
            model_clim[vname][season] = np.mean(season_to_data_for_years[season], axis=0)

    fig = plt.figure()
    nseasons = len(season_to_months)
    gs = GridSpec(3, nseasons + 1, width_ratios=[1.0, ] * nseasons + [0.05, ])
    clevs = model_data_calculator.get_clev_by_name(model_var_name)
    cmap = cm.get_cmap("jet", len(clevs) - 1)
    norm = BoundaryNorm(clevs, len(clevs) - 1)

    clevs_diff = model_data_calculator.get_clev_by_name(model_var_name, delta=True)
    cmap_diff = cm.get_cmap("bwr", len(clevs_diff) - 1)
    norm_diff = BoundaryNorm(clevs_diff, len(clevs_diff) - 1)

    im = None

    #plot model climatology
    col = 0
    for season, data in model_clim[model_var_name].iteritems():
        print col
        ax = fig.add_subplot(gs[0, col])

        if col == 0:
            ax.set_ylabel("Model (ERA40-driven)")

        ax.set_title(season)

        to_plot = np.ma.masked_where(gl_land_sea_mask, data)
        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs[0], vmax=clevs[-1], cmap=cmap,
                                norm=norm)
        basemap.drawmapboundary(fill_color="0.75")
        basemap.drawcoastlines()
        col += 1

    #plot obs climatology
    col = 0
    bias_clim = OrderedDict()
    for season, data in seasonal_obs_data.iteritems():
        ax = fig.add_subplot(gs[1, col])

        if col == 0:
            ax.set_ylabel("Obs (CRU)")

        ax.set_title(season)
        to_plot = cru.interpolate_data_to(data_in=data, lons2d=lons, lats2d=lats, nneighbours=1)

        to_plot = np.ma.masked_where(gl_land_sea_mask, to_plot)
        im = basemap.pcolormesh(x, y, to_plot, vmin=clevs[0], vmax=clevs[-1], cmap=cmap,
                                norm=norm)
        bias_clim[season] = to_plot - model_clim[model_var_name][season]
        basemap.drawmapboundary(fill_color="0.75")
        basemap.drawcoastlines()
        col += 1

    cax = fig.add_subplot(gs[:2, col])
    plt.colorbar(im, ticks=clevs, cax=cax, extend="max" if model_var_name == "PR" else "both")

    col = 0
    for season, data in bias_clim.iteritems():
        ax = fig.add_subplot(gs[2, col])

        if col == 0:
            ax.set_ylabel("Model - Obs")

        im = basemap.pcolormesh(x, y, data, vmin=clevs_diff[0],
                                vmax=clevs_diff[-1], cmap=cmap_diff,
                                norm=norm_diff)
        basemap.drawcoastlines()
        basemap.drawmapboundary(fill_color="0.75")
        col += 1

    plt.colorbar(im, ticks=clevs_diff, extend="both", cax=fig.add_subplot(gs[2, col]))



    img_folder = "../images/validation/cru/{}".format(cru_var_name)
    if not os.path.isdir(img_folder):
        os.makedirs(img_folder)
    img_path = os.path.join(img_folder, "{}-{}.png".format(start_year, end_year))
    fig.savefig(img_path, dpi=100)
    plt.close(fig)



def launch_validation():
    font_size = 25
    fig_size = (25, 18)
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




    #coordinates
    basemap, lons, lats = model_data_calculator.get_basemap_and_coordinates_from_any_file_in(
        folder=ERA40_DRIVEN_DATA_PATH)
    x, y = basemap(lons, lats)
    land_sea_gl_mask = model_data_calculator.get_land_sea_glaciers_mask_from_geophysics_file()

    validate_swe_and_snowdepth(start_year=1999, end_year=2010, basemap=basemap, x=x, y=y,
                               gl_land_sea_mask=land_sea_gl_mask)

    validate_with_cru(start_year=1981, end_year=2010, season_to_months=default_season_name_to_months,
                      basemap=basemap, lons=lons, lats=lats, gl_land_sea_mask=land_sea_gl_mask,
                      model_var_name="TT", cru_var_name="tmp",
                      cru_file="/b2_fs2/huziy/CRUTS3.1/cru_ts_3_10.1901.2009.tmp.dat.nc")

    validate_with_cru(start_year=1981, end_year=2010, season_to_months=default_season_name_to_months,
                      basemap=basemap, lons=lons, lats=lats, gl_land_sea_mask=land_sea_gl_mask,
                      model_var_name="PR", cru_var_name="pre",
                      cru_file="/b2_fs2/huziy/CRUTS3.1/cru_ts_3_10.1901.2009.pre.dat.nc")


if __name__ == "__main__":
    launch_validation()


