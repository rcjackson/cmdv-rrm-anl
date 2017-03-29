import pyart
from matplotlib import pyplot as plt
import yt
import numpy as np
import math
from yt.visualization.volume_rendering.interactive_vr import \
    SceneGraph, BlockCollection, TrackballCamera
from yt.visualization.volume_rendering.interactive_loop import \
    RenderingContext
from yt.visualization.volume_rendering import glfw_inputhook 

grid_data_path = '/home/rjackson/multidop_grids/cf_compliant_grid200601200050.nc'
pyart_grid = pyart.io.read_grid(grid_data_path)
bbox = np.array([[pyart_grid.x['data'][0], pyart_grid.x['data'][-1]],
                 [pyart_grid.y['data'][0], pyart_grid.y['data'][-1]],
                 [pyart_grid.z['data'][0], pyart_grid.z['data'][-1]]])
print(pyart_grid.fields.keys())

units = ['m/s', 'dBZ', 'm/s', 'm/s']
eastward_wind = np.ma.transpose(pyart_grid.fields['eastward_wind']['data'],
                             [2,1,0])
northward_wind = np.ma.transpose(pyart_grid.fields['northward_wind']['data'], 
                              [2,1,0])
upward_air_velocity = np.ma.transpose(pyart_grid.fields['upward_air_velocity']['data'], 
                                   [2,1,0])
reflectivity = np.ma.transpose(pyart_grid.fields['reflectivity']['data'], 
                            [2,1,0])
reflectivity[reflectivity < 0] = np.nan
data = (dict(u = (eastward_wind, 'm/s'),
             v = (northward_wind, 'm/s'),
             w = (upward_air_velocity, 'm/s'),
             dBZ = (reflectivity, 'mm**6/m**3')))
ds = yt.load_uniform_grid(data,
                          eastward_wind.shape,
                          length_unit="m",
                          time_unit="s",
                          velocity_unit="m/s",
                          bbox=bbox,
                          )

sc = yt.create_scene(ds, field=('stream', 'dBZ'), lens_type='perspective')
source = sc[0]
source.set_field("dBZ")
source.set_log(False)
bounds = (0, 60)
camera = sc.add_camera()
#camera.roll(math.pi/2)
#camera.yaw(math.pi/6)
#camera.pitch(math.pi)
#camera.set_focus((0.1,-0.2,0.1))
#camera.set_position((0,0,0))

tf = yt.ColorTransferFunction(bounds)
sc.camera.width = (70000, 'm')
tf.add_layers(5, colormap=pyart.graph.cm.NWSVel)

source.tfh.tf = tf
source.tfh.bounds = bounds
# Start renderer
rc = RenderingContext(1280, 960)

collection = BlockCollection()

dd = ds.all_data()
collection.add_data(dd, "dBZ")

position = (1.0, 1.0, 1.0)
c = TrackballCamera(position=position, focus=ds.domain_center,
                    near_plane=0.1)

callbacks = rc.setup_loop(sc, c)
rl = rc(sc, c, callbacks)

# To make this work from IPython execute:
#
# glfw_inputhook.inputhook_manager.enable_gui("glfw", app=rl)
