import yt
import numpy as np
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, VolumeSource

data=yt.load('disk.out1.00004.athdf')
slc = yt.SlicePlot(data, 'phi', 'Er')
slc.save('image.png')