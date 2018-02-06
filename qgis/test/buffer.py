import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib as plt

#create tree location matrix with values indicating crown radius
A = np.zeros((120,320))
A[60,40] = 1
A[60,80] = 2
A[60,120] = 3
A[60,160] = 4
A[60,200] = 5
A[60,240] = 6
A[60,280] = 7

#plot tree locations
fig = plt.figure()
plt.imshow(A, interpolation='none')
plt.colorbar()

#find unique values
unique_vals = np.unique(A)
unique_vals = unique_vals[unique_vals > 0]

# create circular kernel
def createKernel(radius):
    kernel = np.zeros((2*radius+1, 2*radius+1))
    y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    kernel[mask] = 1
    return kernel

#apply binary dilation sequentially to each unique crown radius value
C = np.zeros(A.shape).astype(bool)
for k, radius in enumerate(unique_vals):
    B = ndimage.morphology.binary_dilation(A == unique_vals[k], structure=createKernel(radius))
    C = C | B #combine masks

#plot resulting mask
fig = plt.figure()
plt.imshow(C, interpolation='none')
plt.show()