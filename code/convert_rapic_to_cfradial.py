import os

out_dir = '/home/rjackson/data/radar/berr/2009'
in_dir = '/lcrc/group/earthscience/radar/stage/radar_disk_two/berr_rapic/berr_2009/'

for root, dirs, files in os.walk(in_dir, topdown=False):
    for name in dirs:
        os.chdir(os.path.join(root, name))
        print('In directory ' + os.path.join(root, name))
        os.system('RadxConvert -f *.rapic -v -outdir ' + out_dir)
  



