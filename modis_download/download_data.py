from multiprocessing import Pool
import os

__author__ = 'huziy'

"""
The script is used to selectively download mmodis albedo datta for selected years.
Data documentation: https://lpdaac.usgs.gov/products/modis_products_table/mcd43c3
Here I am interested in total surface albedo
"""

base_url = "http://e4ftl01.cr.usgs.gov/MOTA/MCD43C3.005/"

import requests
import re


def get_file_links(remote_folder):
    response = requests.get(remote_folder)
    groups = re.findall(r"href=\"(?P<filename>\w+\.\w+\.\w+.\w+\.hdf)\"", response.content)
    if len(groups) != 1:
        print "{} has more than {} hdf file inside ...".format(remote_folder, len(groups))
    return [os.path.join(remote_folder, group) for group in groups]


def download_for_day(args):
    """
    Download data for a day
    :param args:
    :return:
    """
    folder_name, target_dir = args
    day_dir_local = os.path.join(target_dir, folder_name)
    if not os.path.isdir(day_dir_local):
        os.makedirs(day_dir_local)

    day_dir_remote = os.path.join(base_url, folder_name)
    for file_link in get_file_links(day_dir_remote):
        fname = os.path.basename(file_link)
        fpath_local = os.path.join(day_dir_local, fname)


        if os.path.isfile(fpath_local):
            remote_size = int(requests.head(file_link).headers["content-length"])
            local_size = os.path.getsize(fpath_local)

            print "local-size={}; remote-size={}".format(local_size, remote_size)
            if local_size >= remote_size:
                return

        print "downloading {}".format(file_link)
        with open(fpath_local, "wb") as f:
            response = requests.get(file_link)
            for chunk in response.iter_content():
                f.write(chunk)




def main(start_year=2001, end_year=2010,
         target_dir="/b2_fs2/huziy/PythonProjects/DATA/modis_surface_albedo"):

    response = requests.get(base_url)
    assert isinstance(response, requests.Response)
    folder_names = re.findall(r"\d+\.\d+\.\d+", response.content)

    years = [int(fname.split(".")[0]) for fname in folder_names]
    folder_names = [fname for fname, year in zip(folder_names, years) if start_year <= year <= end_year]
    print len(folder_names)
    print folder_names[:10]

    pool = Pool(processes=10)

    target_dirs = [target_dir, ] * len(folder_names)
    pool.map(download_for_day, zip(folder_names, target_dirs))


if __name__ == "__main__":
    main()