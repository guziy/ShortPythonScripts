from datetime import datetime

__author__ = 'huziy'


def get_year_and_month_from_name(fname):
    try:
        d = datetime.strptime(fname.split("_")[-1], "%Y%m")
        return d.year, d.month
    except ValueError:
        return -1, -1
