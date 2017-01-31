import glob

in_dir = '/home/rjackson/data/radar/'
out_dir = '/home/rjackson/data/radar/berr/'

def split_berrima_file(file_name):
    # Look for /IMAGE and /IMAGEEND in file
    rapic_file = open(file_name, 'rb')
    string = b'!'
    while string:
        while not b'/IMAGE' in string:
            string = rapic_file.readline()
            if string == b'':
                break
            print(string)
        if string == b'':
            break
        # Get start time and write to destination file
        date_string = string[13:26]
        year_string = date_string[0:2]
        month_string = date_string[2:4]
        day_string = date_string[4:6]
        hour_string = date_string[6:8]
        minute_string = date_string[8:10]
        print('Found scan at ' + 
              str(year_string, 'utf-8') +
              '-' +
              str(month_string, 'utf-8') +
              '-' +
              str(day_string, 'utf-8') + 
              ' ' +
              str(hour_string, 'utf-8') + 
              ':' +
              str(minute_string, 'utf-8'))
        # Output to split file
        out_file_name = (out_dir + 
                        'BerrimaVol' + 
                         str(year_string, 'utf-8') + 
                         str(month_string, 'utf-8') + 
                         str(day_string, 'utf-8') + 
                         str(hour_string, 'utf-8') + 
                         str(minute_string, 'utf-8') + 
                        '.rapic')
        out_rapic_file = open(out_file_name, 'wb+')
        out_rapic_file.write(string)
        while not b'/IMAGEEND' in string:
            string = rapic_file.readline()          
            if string == b'':
                break
            out_rapic_file.write(string)
  
        out_rapic_file.close()
        if string == b'':
            break
        string = rapic_file.readline()

# Loop over each file in directory
file_list = glob.glob(in_dir + '*.rapic')
print(file_list)
for file_name in file_list:
    print('Processing ')
    split_berrima_file(file_name)
