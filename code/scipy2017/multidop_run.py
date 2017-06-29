from multidop_time import do_multidop_for_time
from distributed import Client, LocalCluster, wait
import time_procedures
from time import sleep
import time
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import sys
import os

""" Get the radars from given time.
Input the range of dates and time wanted for the collection of images """
start_year = int(sys.argv[1])
start_month = int(sys.argv[2])
start_day = int(sys.argv[3])
end_year = int(sys.argv[4])
end_month = int(sys.argv[5])
end_day = int(sys.argv[6])

# serial = 0 run in parallel, = 1 run in serial
serial = 0
times = time_procedures.get_radar_times_cpol(start_year, start_month,
                                             start_day, 19, 0, end_year,
                                             end_month, end_day, 0, 2)
if(serial == 0):
    # Initalize the cluster. Adjust the number of workers to your liking.
    Cluster = LocalCluster(n_workers=4, processes=False)
    client = Client(Cluster)
    # Map the calls to multidop onto the workers
    the_futures = client.map(do_multidop_for_time, times[0])
    wait(the_futures)
else:
    for timer in times[0]:
        do_multidop_for_time(timer)
