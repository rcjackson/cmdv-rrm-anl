import glob

data_file_path = '/lcrc/group/earthscience/radar/stage/disk_one/cpol_lassen/cpol_0405/'

out_file = open('list_of_files0405', 'w+')
i = 0

files = glob.glob(data_file_path + '*PPI.lassen')

for lines in files:
    out_file.write(lines + '\n')
        
out_file.close()




