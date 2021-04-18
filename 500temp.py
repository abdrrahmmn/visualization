from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature
from cartopy.io.shapereader import Reader
from matplotlib.offsetbox import  OffsetImage
import matplotlib.image as image
import datetime

from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)

# Open the NetCDF file
ncfile = Dataset("/home/kodlama/data/wrfout_d01_2020-11-25.nc")

# shapefile oku, il sınırları

fname = "/home/kodlama/visualization/shapefiles/gadm36_TUR_1.shp"

# Extract the pressure, geopotential height, and wind variables
# Burada istediğimiz değişkenleri çekebiliriz, tablosu var.
timeg=1

p = getvar(ncfile, "pressure",timeidx=timeg)
slp = getvar(ncfile, "slp",units='hPa',timeidx=timeg)
z=getvar(ncfile,"z",units="m",timeidx=timeg)
t= getvar(ncfile,'tc',timeidx=timeg)

time=g_times.get_times(ncfile,timeidx=timeg,meta=False)
time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')

ht500=vinterp(ncfile,
               field=z,
               vert_coord="p",
               interp_levels=[500],
               extrapolate=True,
               field_type="z",
               log_p=True)
t500=vinterp(ncfile,
               field=t,
               vert_coord="p",
               interp_levels=[500],
               extrapolate=True,
               field_type="tc",
               log_p=True)

smooth_slp = smooth2d(slp, 5, cenweight=4)
smooth_t500 = smooth2d(t500, 3, cenweight=4)

# Get the lat/lon coordinatesç
# Verinin enlem ve boylamini otomatik almasi icin yazildi.
lats, lons = latlon_coords(slp)

# Get the map projection information, verinin projeksiyon methodunu belirlemek icin yazildi.
cart_proj = get_cartopy(slp)

# Create the figure, olusturulacak goruntunun boyut ayari.
fig = plt.figure(figsize=(10,6))

ax = plt.axes(projection=cart_proj)
logo=plt.imread("logo.png")

ax.figure.figimage(logo, 380, 400,alpha=0.7)
# Download and add the states and coastlines, icindeki parametreler degisince harita cizimi degisir.
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_0_boundary_lines_land")
# Yukarida belirlenen ayarlari figure eklemek icin add_feature kullanildi.
ax.add_feature(states, linewidth=0.4, edgecolor="dimgray")
# Sahil kenarlarini cizmeye yarar.
ax.add_feature(COASTLINE,linewidth=0.5,edgecolor="dimgray")
shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.2)
ax.add_feature(shape_feature)

# Cizilecek contour lari plt.contour ile yapiyoruz. Mesela to_np(lons) numpy arrayine ceviriyor boylam degerlerini. Bunlarin sirasinin onemi var mi bilmiyorum. transform onemli, cartopy kullandik, basemap de kullanilabilirdi.
contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                        colors="black",
                        transform=crs.PlateCarree(),linewidths=1)
# Burada contour label ayarlari yapildi.
plt.clabel(contours, inline=1, fontsize=10, fmt="%i",colors="black")
levels = np.arange(4760,6040,40)
ht500contour= plt.contourf(to_np(lons), to_np(lats),to_np(ht500[0]),extend="both",cmap="rainbow",
            transform=crs.PlateCarree(),levels=levels)
cbar= plt.colorbar(ht500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)

cbar.set_label("500 hPa Jeopotansiyel Yükseklik (metre)",labelpad=-40)

tcontour=plt.contour(to_np(lons), to_np(lats), to_np(smooth_t500[0]),colors="black",
                        transform=crs.PlateCarree(),linewidths=0.5,linestyles="--")
plt.clabel(tcontour, inline=1, fontsize=8, fmt="%i",colors="black")

# Set the map bounds, figurde eksenlerin sinirlari ve degerlerini alabilmek icin yazildi.
ax.set_xlim(cartopy_xlim(slp))
ax.set_ylim(cartopy_ylim(slp))

#title=ax.set_title('500 hPa Sıcaklığı (°C), 500 hPa Jeopotansiyel Yüksekliği (m), Deniz Seviyesi Basıncı (mb)',backgroundcolor='black',pad=10)
plt.title("500 hPa Sıcaklığı (°C)\n500 hPa Jeopotansiyel Yüksekliği (m)\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
# Figuru kaydetmek icin plt.savefig kullanildi.
plt.savefig('500hpatmpgph.png',dpi=300,transparent=False)
