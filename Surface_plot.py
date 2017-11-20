from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pylab as pl
from PIL import Image
from PIL import ImageChops
from PIL import ImageFilter
import numpy as np
import pylab
import scipy.misc
img1 = Image.open('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160518/WiTEc/test13.bmp').convert('L')
fil1 = img1.filter(ImageFilter.MedianFilter(7))
img2 = Image.open('/home/dwi/ownCloud/Anlabshared/Dwi/Data and Analysis/160518/WiTEc/test14.bmp').convert('L')
fil2 = img2.filter(ImageFilter.MedianFilter(7))
img = ImageChops.subtract(fil2, fil1, 1) 
z   = np.asarray(img)
mydata1 = z[::1,::1]
mydata = scipy.misc.imresize(mydata1,100, interp='cubic')
fig = pl.figure()
#ax1 = fig.add_subplot(1,2,1)
#pl.imshow(mydata1)
#ax1.set_title('2D')
#print z
ax2 = fig.gca(projection='3d')
x,y = np.mgrid[:mydata.shape[0],:mydata1.shape[1]]
ax2.plot_surface(x,y,mydata,cmap=pl.cm.coolwarm,rstride=1,cstride=1,linewidth=0.,antialiased=False)
ax2.set_title('3D')
ax2.set_zlim3d(0,150)
#pl.colorbar()
#pl.colorbar(shrink=0.5, aspect=5)
pl.show()
