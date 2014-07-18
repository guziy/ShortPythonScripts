__author__ = 'huziy'



#The data is on the CMG (climate modeling grid)
#Grid description can be found here: http://modis-land.gsfc.nasa.gov/MODLAND_grid.html

from pyresample.geometry import AreaDefinition

dx = dy = 0.05
nx = int(360.0 / dx)
ny = int(180.0 / dy)

data_shape = (ny, nx)
area_extent = (-180.0, 90.0, 180.0, -90.0)

proj_dict = {
    "proj": "latlong",
    "ellps": "WGS84",
    "datum": "WGS84"
}

pr_area_definition = AreaDefinition("mcd43c3", "mcd43c3", "latlong", proj_dict, nx, ny, area_extent=area_extent)

#print nx, ny, pr_area_definition

