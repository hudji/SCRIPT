from landlab import RasterModelGrid
import numpy as np
grid = RasterModelGrid((4, 5))#
elevation = np.array([0., 0., 0., 0., 0, 1., 1., 1., 1., 1, 2., 2., 2., 2., 2, 3., 3., 3., 3., 3])