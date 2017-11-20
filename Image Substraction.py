import numpy as np
import matplotlib.pyplot as pl
import matplotlib.image as im
from mpl_toolkits.mplot3d import Axes3D
import scipy.misc  
fig = pl.figure()
#ax = fig.gca(projection='3d')
image1 = im.imread ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160518/WiTEc/test13.bmp')
image2 = im.imread ('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160518/WiTEc/test14.bmp')
#print image1
imsubs = image2 - image1
imres = scipy.misc.imresize (imsubs, 30, interp = 'cubic')
#Yaxis = np.array[len(imsubs)]
#Xaxis = np.array[len(imsubs[0])]
#X, Y = np.mgrid [0:mydata.shape[0], 0:mydata.shape[1]]
#Z = imarr.reshape(Y.size,X.size)
#ax.plot_surface(X, Y, mydata, rstride=1, cstride=1, cmap=pl.cm.gray, linewidth=0)
pl.imshow (imres, interpolation = "bicubic")
#print imsubs
#print shape
#print Y
#xx = np.len[imres]
#print xx
pl.colorbar()
pl.show()