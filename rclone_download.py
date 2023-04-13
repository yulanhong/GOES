import subprocess
import os
from datetime import datetime, date
import numpy as np

year = '2022'

month='02'

days = np.arange(4)+1

for day in days:

    savedir = "/data/sfrdata/yulanh/data/GOES/"+year+month+str(day).zfill(2)
    print(savedir)
    if os.path.exists(savedir) == False:
        os.system('mkdir '+savedir)

    for hour in np.arange(24):

        hh=str(hour).zfill(2)
        print(hh)

        day_of_year = date(int(year),int(month),day).timetuple().tm_yday

        downdir = "mypubAWS:noaa-goes16/ABI-L1b-RadF/"+year+\
                "/"+str(day_of_year).zfill(3)+"/"+hh+"/"

        files = subprocess.check_output('rclone ls '+downdir, shell=True)

        files = files.decode()

        files = files.split('\n')
        files.remove('')
        files =  [i.split(" ")[-1] for i in files]
        #print(files)
# Print the list of directories
        os.system("rclone copy "+downdir+'  '+savedir)
