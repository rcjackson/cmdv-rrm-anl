from multidop_time import do_multidop_for_time
from ipyparallel import Client
import time_procedures
from time import sleep
import time
import numpy as np

# Get the radars from given time.
# Input the range of dates and time wanted for the collection of images
start_year = 2006
start_day = 19
start_month = 1
start_hour = 23
start_minute = 00
start_second = 0

end_year = 2006
end_month = 1
end_day = 20
end_hour = 1
end_minute = 0
end_second = 0


# Make the sounding file for the input
times = time_procedures.get_radar_times_cpol(start_year, start_month, 
                                             start_day, start_hour, 
                                             start_minute, end_year,
                                             end_month, end_day, end_hour, 
                                             end_minute, minute_interval=0)


# Get iPython cluster
#state = 0
#while state == 0:
#    try:
#        My_Cluster = Client()
#        My_View = My_Cluster[:]
#        state = 1
#    except:
#        state = 0
#        print('Cluster not ready for me')
#        sleep(10)

#Turn off blocking so all engines can work async
#My_View.block = False

#on all engines do an import of Py-ART
#My_View.execute('import matplotlib')
#My_View.execute('matplotlib.use("agg")')
#My_View.execute('import pyart')
#My_View.execute('import numpy as np')
t1 = time.time()


#for timer in times:
#   do_multidop_for_time(timer)

#Map the code and input to all workers
#result = My_View.map_async(do_multidop_for_time, times)


#Reduce the result to get a list of output
qvps = result.get()
tt = time.time() - t1
print(tt)
print(tt/len(times))


    
