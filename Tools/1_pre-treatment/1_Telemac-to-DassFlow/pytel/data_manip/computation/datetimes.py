#!/usr/bin/env python3
""" Datetimes computations """
from datetime import datetime, timedelta
import numpy as np


def compute_datetimes(times,
                      initial_date='2019/01/01 00:00:00',
                      time_format='%Y/%m/%d %H:%M:%S'):
    """
    Convert times from seconds to datetimes

    @param times (list) List of time in seconds
    @param initial_date (datetime, str or np.array) initial date
    (default:'01/01/2019 00:00:00')
    @param time_format (str) time format (default: '%Y/%m/%d %H:%M:%S')
    """
    if isinstance(initial_date, str):
        ini_datetime_obj = datetime.strptime(initial_date, time_format)
    elif isinstance(initial_date, np.ndarray):
        datetime_string = '%i/%i/%i %i:%i:%i' % (
            initial_date[0], initial_date[1], initial_date[2],
            initial_date[3], initial_date[4], initial_date[5])
        ini_datetime_obj = datetime.strptime(datetime_string, time_format)
    else:
        ini_datetime_obj = initial_date

    ntimes = len(times)
    datetimes = [ini_datetime_obj for t in range(ntimes)]
    for time in range(ntimes):
        datetimes[time] += timedelta(seconds=times[time])

    return datetimes
