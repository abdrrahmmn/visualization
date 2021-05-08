import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import salem, datetime
from matplotlib.offsetbox import  OffsetImage
import matplotlib.image as image
from cmapcustom import colorgetcustom
# ncfiledir=input('Ncfile addres')
# fnamedir=input('fname')

customcmap=colorgetcustom()

ncfiledir = "/home/durmaz/wrf/Build_WRF/WRF/run/wrfout_d01_2021-05-02_00:00:00"
fnamedir = "/home/durmaz/visualization/shapefiles/gadm36_TUR_1.shp"

def vor_500h(ncfiledir,fnamedir,timeg):
    from wrf import smooth2d

    ds = salem.open_wrf_dataset(ncfiledir).isel(time=timeg)
    timestart = ds.xtime
    timestart = np.datetime_as_string(timestart.values)
    timestart = datetime.datetime.strptime(str(timestart),'%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=timeg)

    timeout = ds.time
    timeout = np.datetime_as_string(timeout.values)
    timeout = datetime.datetime.strptime(str(timeout),'%Y-%m-%dT%H:%M:%S.%f000')

    # Vorticty 500 height mb
    from netCDF4 import Dataset
    from wrf import getvar,vinterp
    ncfile = Dataset(str(ncfiledir))
    vor = getvar(ncfile,"avo",timeidx=timeg)
    ht500 = ds.salem.wrf_plevel('Z', levels=500)
    vor500=vinterp(ncfile,
                field=vor,
                vert_coord="p",
                interp_levels=[500],
                extrapolate=True,
                log_p=True)
    smooth_ht500 = smooth2d(ht500, 2, cenweight=4)

    fig, ax = plt.subplots(figsize=(10,6))
    # Read shape file


    # plot the salem map background, make countries in grey
    smap = ds.salem.get_map(countries=False)
    smap.set_shapefile(countries=True, color='black')
    smap.set_lonlat_contours(interval=0)
    # Add shapefile
    smap.set_shapefile(fnamedir,zorder=2)
    smap.plot(ax=ax)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    logo=plt.imread("logo.png")
    plt.figimage(logo, 380, 400,alpha=0.7)

    xx, yy = smap.grid.transform(ht500.west_east.values, ht500.south_north.values,
                                 crs=ht500.salem.grid.proj)
    contours = ax.contour(xx, yy, smooth_ht500/10, colors="dimgray", linewidths=1,levels=np.arange(476,604,4))
    ax.clabel(contours, inline=1, fontsize=10, fmt="%i",colors="white")

    levels = np.arange(-12,13,2)
    vor500contour= ax.contourf(xx, yy, vor500[0],extend="both",cmap="coolwarm",
                    levels=levels)
    cbar= fig.colorbar(vor500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("500 hPa Vortisiti ($10^{-5}$.s-1)",labelpad=-40)
    plt.title("500 hPa Vortisiti ($10^{-5}$.$s^{-1}$) \n500 hPa Jeopotansiyel Yüksekliği (dam)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(timestart,timeout),loc="right",pad=10,fontsize=10)
    plt.savefig('vor_500h_{}.png'.format(timeout),dpi=300,transparent=False)

def pre_slp_wind(ncfiledir,fnamedir):
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature,COLORS,LAND,OCEAN,BORDERS
    from cartopy.io.shapereader import Reader
    from netCDF4 import Dataset
    from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    alltimes= getvar(ncfile, "times",ALL_TIMES)

    for timeindex in np.arange(3,len(alltimes),3):

        slp = getvar(ncfile, "slp",timeidx=timeindex,units='hPa')
        RAINNC = getvar(ncfile, "RAINNC",timeidx=timeindex)
        RAINC = getvar(ncfile, "RAINC",timeidx=timeindex)
        RAINSH = getvar(ncfile, "RAINSH",timeidx=timeindex)
        SNOWNC = getvar(ncfile, "SNOWNC",timeidx=timeindex)
        GRAUPELNC = getvar(ncfile, "GRAUPELNC",timeidx=timeindex)
        HAILNC = getvar(ncfile, "HAILNC",timeidx=timeindex)
        u,v= getvar(ncfile,'uvmet10',timeidx=timeindex)

        pre_total = to_np(RAINC[:]) + to_np(RAINNC[:]) + to_np(SNOWNC[:]) + to_np(GRAUPELNC[:]) + to_np(HAILNC[:]+to_np(RAINSH[:]))

        RAINNC = getvar(ncfile, "RAINNC",timeidx=timeindex-3)
        RAINC = getvar(ncfile, "RAINC",timeidx=timeindex-3)
        RAINSH = getvar(ncfile, "RAINSH",timeidx=timeindex-3)
        SNOWNC = getvar(ncfile, "SNOWNC",timeidx=timeindex-3)
        GRAUPELNC = getvar(ncfile, "GRAUPELNC",timeidx=timeindex-3)
        HAILNC = getvar(ncfile, "HAILNC",timeidx=timeindex-3)

        ex_pre_total = to_np(RAINC[:]) + to_np(RAINNC[:]) + to_np(SNOWNC[:]) + to_np(GRAUPELNC[:]) + to_np(HAILNC[:]+to_np(RAINSH[:]))

        pre_total=pre_total-ex_pre_total
        smooth_slp = smooth2d(slp[:], 50, cenweight=4)

        time=g_times.get_times(ncfile,timeidx=timeindex,meta=False)
        time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

        baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
        baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')


        lats, lons = latlon_coords(slp)

        cart_proj = get_cartopy(slp)

        fig = plt.figure(figsize=(10,6))

        ax = plt.axes(projection=cart_proj)
        logo=plt.imread("logo.png")

        ax.figure.figimage(logo, 380, 400,alpha=0.7)
        
        ax.add_feature(LAND, linewidth=0, edgecolor="dimgray",zorder=1)
        ax.add_feature(OCEAN,linewidth=0,edgecolor="dimgray",zorder=1)
        ax.add_feature(BORDERS,linewidth=0.4,edgecolor="dimgray",zorder=2)
        ax.add_feature(COASTLINE,linewidth=0.8,edgecolor="dimgray",zorder=2)

        shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                        crs.PlateCarree(), facecolor='none',linewidth=0,zorder=1)
        ax.add_feature(shape_feature)

        shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                        crs.PlateCarree(),facecolor='none',edgecolor="dimgray",linewidth=0.2,zorder=2)
        ax.add_feature(shape_feature)

        levels = np.arange(0,61,3)
        precip= plt.contourf(to_np(lons), to_np(lats),np.ma.masked_less((pre_total),0.2),extend="both",cmap=customcmap,
                    transform=crs.PlateCarree(),levels=levels,zorder=2)
        cbar= plt.colorbar(precip,fraction=0.071,orientation="horizontal",pad=.07,aspect=50,format="%i",ticks=np.arange(0,61,5))
        cbar.set_label("Yağış (mm) ",labelpad=-40)

        plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
        to_np(u[::10, ::10]), to_np(v[::10, ::10]),
        transform=crs.PlateCarree(), length=4,zorder=2,linewidth=0.4)

        slpcontours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                    colors="maroon", transform=crs.PlateCarree(),linewidths=1,zorder=2)
        plt.clabel(slpcontours, inline=1, fontsize=5, fmt="%i",colors="maroon")

        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))

        plt.title("3 Saatlik Yağış (mm)\n10m Rüzgar (kt)\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
        plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

        plt.savefig('pre_slp_wind_{}.png'.format(time),dpi=300,transparent=False)

def pre(ncfiledir,fnamedir):
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature,COLORS,LAND,OCEAN,BORDERS
    from cartopy.io.shapereader import Reader
    from netCDF4 import Dataset
    from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    alltimes= getvar(ncfile, "times",ALL_TIMES)

    for timeindex in np.arange(3,len(alltimes),3):

        slp = getvar(ncfile, "slp",timeidx=timeindex,units='hPa')
        RAINNC = getvar(ncfile, "RAINNC",timeidx=timeindex)
        RAINC = getvar(ncfile, "RAINC",timeidx=timeindex)
        RAINSH = getvar(ncfile, "RAINSH",timeidx=timeindex)
        SNOWNC = getvar(ncfile, "SNOWNC",timeidx=timeindex)
        GRAUPELNC = getvar(ncfile, "GRAUPELNC",timeidx=timeindex)
        HAILNC = getvar(ncfile, "HAILNC",timeidx=timeindex)
        

        pre_total = to_np(RAINC[:]) + to_np(RAINNC[:]) + to_np(SNOWNC[:]) + to_np(GRAUPELNC[:]) + to_np(HAILNC[:]+to_np(RAINSH[:]))

        RAINNC = getvar(ncfile, "RAINNC",timeidx=timeindex-3)
        RAINC = getvar(ncfile, "RAINC",timeidx=timeindex-3)
        RAINSH = getvar(ncfile, "RAINSH",timeidx=timeindex-3)
        SNOWNC = getvar(ncfile, "SNOWNC",timeidx=timeindex-3)
        GRAUPELNC = getvar(ncfile, "GRAUPELNC",timeidx=timeindex-3)
        HAILNC = getvar(ncfile, "HAILNC",timeidx=timeindex-3)

        ex_pre_total = to_np(RAINC[:]) + to_np(RAINNC[:]) + to_np(SNOWNC[:]) + to_np(GRAUPELNC[:]) + to_np(HAILNC[:]+to_np(RAINSH[:]))

        pre_total=pre_total-ex_pre_total
        smooth_slp = smooth2d(slp[:], 50, cenweight=4)

        time=g_times.get_times(ncfile,timeidx=timeindex,meta=False)
        time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

        baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
        baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')


        lats, lons = latlon_coords(slp)

        cart_proj = get_cartopy(slp)

        fig = plt.figure(figsize=(10,6))

        ax = plt.axes(projection=cart_proj)
        logo=plt.imread("logo.png")

        ax.figure.figimage(logo, 380, 400,alpha=0.7)
        ax.add_feature(LAND, linewidth=0, edgecolor="dimgray",zorder=1)
        ax.add_feature(OCEAN,linewidth=0,edgecolor="dimgray",zorder=1)
        ax.add_feature(BORDERS,linewidth=0.4,edgecolor="dimgray",zorder=2)
        ax.add_feature(COASTLINE,linewidth=0.8,edgecolor="dimgray",zorder=2)

        shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                        crs.PlateCarree(), facecolor='none',linewidth=0,zorder=1)
        ax.add_feature(shape_feature)

        shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                        crs.PlateCarree(),facecolor='none',edgecolor="dimgray",linewidth=0.2,zorder=2)
        ax.add_feature(shape_feature)

        levels = np.arange(0,61,3)
        precip= plt.contourf(to_np(lons), to_np(lats),np.ma.masked_less((pre_total),0.2),extend="both",cmap=customcmap,
                    transform=crs.PlateCarree(),levels=levels,zorder=2)
        cbar= plt.colorbar(precip,fraction=0.071,orientation="horizontal",pad=.07,aspect=50,format="%i",ticks=np.arange(0,61,5))
        cbar.set_label("Yağış (mm) ",labelpad=-40)


        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        plt.title("3 Saatlik Yağış (mm)",loc="left",pad=10,fontsize=12)
        plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

        plt.savefig('pre_{}.png'.format(time),dpi=300,transparent=False)

def slp_500h(ncfiledir,fnamedir,timeg):
    from wrf import smooth2d

    ds = salem.open_wrf_dataset(ncfiledir).isel(time=timeg)
    timestart = ds.xtime
    timestart = np.datetime_as_string(timestart.values)
    timestart = datetime.datetime.strptime(str(timestart),'%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=timeg)

    timeout = ds.time
    timeout = np.datetime_as_string(timeout.values)
    timeout = datetime.datetime.strptime(str(timeout),'%Y-%m-%dT%H:%M:%S.%f000')
    slp = ds.SLP
    ht500 = ds.salem.wrf_plevel('Z', 500)
    smooth_slp = smooth2d(slp,100, cenweight=4)

    # Bunu subplots ne öğren ?
    fig, ax = plt.subplots(figsize=(10,6))



    # plot the salem map background, make countries in grey
    smap = ds.salem.get_map(countries=False)
    smap.set_shapefile(countries=True, color='grey')
    smap.set_lonlat_contours(interval=0)
    # Add shapefile
    smap.set_shapefile(fnamedir,zorder=2)
    smap.plot(ax=ax)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)


    logo=plt.imread("logo.png")
    plt.figimage(logo, 380, 400,alpha=0.7)

    xx, yy = smap.grid.transform(slp.west_east.values, slp.south_north.values,
                                 crs=slp.salem.grid.proj)
    levels = np.arange(950,6040,40)
    contours = ax.contour(xx, yy, smooth_slp,colors="black",linewidths=1,levels=np.arange(950,1100,2))
    ax.clabel(contours, inline=1, fontsize=7, fmt=lambda p: f'{str(int(p))[-2:]}')

    levels = np.arange(476,604,4)
    ht500contour= ax.contourf(xx, yy,ht500/10,extend="both",cmap="rainbow",levels=levels)
    cbar= fig.colorbar(ht500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("500 hPa Jeopotansiyel Yükseklik (dam)",labelpad=-40)

    plt.title("500 hPa Jeopotansiyel Yüksekliği (m)\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
    plt.title(f'Başlangıç: {timestart} UTC\nGeçerli: {timeout} UTC',loc="right",pad=10,fontsize=10)

    fig.savefig('slp_500h_{}.png'.format(timeout),dpi=300,transparent=False) 

def slp_wind_t2(ncfiledir,fnamedir,timeg):
    from wrf import smooth2d

    ds = salem.open_wrf_dataset(ncfiledir).isel(time=timeg)
    timestart = ds.xtime
    timestart = np.datetime_as_string(timestart.values)
    timestart = datetime.datetime.strptime(str(timestart),'%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=timeg)

    timeout = ds.time
    timeout = np.datetime_as_string(timeout.values)
    timeout = datetime.datetime.strptime(str(timeout),'%Y-%m-%dT%H:%M:%S.%f000')
    #Temperature at 2m and wind slp
    t2c = ds.T2C
    u10 = ds.U10 
    v10 = ds.V10
    slp = ds.SLP
    smooth_slp = smooth2d(slp, 100, cenweight=4)
    # Bunu subplots ne öğren ?
    fig, ax = plt.subplots(figsize=(10,6))
    # Read shape file

    # plot the salem map background, make countries in grey
    smap = ds.salem.get_map(countries=False)
    smap.set_shapefile(countries=True, color='grey')
    smap.set_lonlat_contours(interval=0)
    # Add shapefile
    smap.set_shapefile(fnamedir,zorder=2)
    smap.plot(ax=ax)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    logo=plt.imread("logo.png")
    plt.figimage(logo, 380, 400,alpha=0.7)

    xx, yy = smap.grid.transform(t2c.west_east.values, t2c.south_north.values,
                                 crs=t2c.salem.grid.proj)

    levels = np.arange(950,1100,2)
    slpcontours = ax.contour(xx,yy,smooth_slp,
                                colors="maroon",linewidths=1,levels=levels)
    ax.clabel(slpcontours, inline=1, fontsize=10, fmt=lambda p: f'{str(int(p))[-2:]}',colors="maroon")
    levels = np.arange(-30,44,2)
    t2mcontour= ax.contourf(xx, yy, t2c, extend="both",cmap="jet",
                    levels=levels)

    cbar= fig.colorbar(t2mcontour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("2 metre Sıcaklığı (°C) ",labelpad=-40)
    ax.barbs(xx[::12], yy[::12],
                u10[::12, ::12], v10[::12, ::12],
                 length=5)

    plt.title("Deniz Seviyesi Basıncı (hPa)\n2 metre Sıcaklığı (°C)\n10 metre rüzgar (kt)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(timestart,timeout),loc="right",pad=10,fontsize=10)

    fig.savefig('slp_wind_t2_{}.png'.format(timeout),dpi=300,transparent=False)

def wind_t2(ncfiledir,fnamedir,timeg):
    ds = salem.open_wrf_dataset(ncfiledir).isel(time=timeg)
    timestart = ds.xtime
    timestart = np.datetime_as_string(timestart.values)
    timestart = datetime.datetime.strptime(str(timestart),'%Y-%m-%dT%H:%M:%S.%f000') - datetime.timedelta(hours=timeg)

    timeout = ds.time
    timeout = np.datetime_as_string(timeout.values)
    timeout = datetime.datetime.strptime(str(timeout),'%Y-%m-%dT%H:%M:%S.%f000')
    # Temperature at 2m and wind
    t2c = ds.T2C
    u10 = ds.U10 
    v10 = ds.V10 
    # Bunu subplots ne öğren ?
    fig, ax = plt.subplots(figsize=(10,6),dpi=300)
    # Read shape file


    # plot the salem map background, make countries in grey
    smap = ds.salem.get_map(countries=False)
    smap.set_shapefile(countries=True, color='black')
    smap.set_lonlat_contours(interval=0)
    # Add shapefile
    smap.set_shapefile(fnamedir,zorder=2)
    smap.plot(ax=ax)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    logo=plt.imread("logo.png")
    plt.figimage(logo, 380, 400,alpha=0.7)

    xx, yy = smap.grid.transform(t2c.west_east.values, t2c.south_north.values,
                                 crs=t2c.salem.grid.proj)
    levels = np.arange(-30,44,2)
    t2mcontour= ax.contourf(xx, yy, t2c, extend="both",cmap="jet",
                    levels=levels)

    cbar= fig.colorbar(t2mcontour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("2 metre Sıcaklığı (°C) ",labelpad=-40)
    ax.barbs(xx[::12], yy[::12],
                u10[::12, ::12], v10[::12, ::12],
                 length=5)
    plt.title("2 metre Sıcaklığı (°C)\n10 metre rüzgar (kt)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(timestart,timeout),loc="right",pad=10,fontsize=10)
    fig.savefig('wind_t2_{}.png'.format(timeout),dpi=300,transparent=False)   

def rh850_wind850(ncfiledir,fnamedir,timeg):
    from netCDF4 import Dataset
    from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature
    from cartopy.io.shapereader import Reader
    from matplotlib.offsetbox import  OffsetImage
    import matplotlib.image as image

    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)

    time=g_times.get_times(ncfile,timeidx=timeg,meta=False)
    time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

    baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
    baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')

    # Extract the pressure, geopotential height, and wind variables
    # Burada istediğimiz değişkenleri çekebiliriz, tablosu var.
    rh = getvar(ncfile,"rh",timeidx=timeg)
    u,v= getvar(ncfile,"uvmet",timeidx=timeg)

    rh850=vinterp(ncfile,
                field=rh,
                vert_coord="p",
                interp_levels=[850],
                extrapolate=True,
                log_p=True)

    u850=vinterp(ncfile,
                field=u,
                vert_coord="p",
                interp_levels=[850],
                extrapolate=True,
                log_p=True)
    v850=vinterp(ncfile,
                field=v,
                vert_coord="p",
                interp_levels=[850],
                extrapolate=True,
                log_p=True)



    # Get the lat/lon coordinatesç
    # Verinin enlem ve boylamini otomatik almasi icin yazildi.
    lats, lons = latlon_coords(rh850)

    # Get the map projection information, verinin projeksiyon methodunu belirlemek icin yazildi.
    cart_proj = get_cartopy(rh850)

    # Create the figure, olusturulacak goruntunun boyut ayari.
    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)
    # Download and add the states and coastlines, icindeki parametreler degisince harita cizimi degisir.
    states = NaturalEarthFeature(category="cultural", scale="10m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")
    # Yukarida belirlenen ayarlari figure eklemek icin add_feature kullanildi.
    ax.add_feature(states, linewidth=0.8, edgecolor="dimgray")
    # Sahil kenarlarini cizmeye yarar.
    ax.add_feature(COASTLINE,linewidth=0.9,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.5)
    ax.add_feature(shape_feature)



    rh850contour=plt.contourf(to_np(lons), to_np(lats), to_np(rh850[0]),transform=crs.PlateCarree(),cmap="Purples",levels=np.arange(0,101,10),extend="both")
    cbar= plt.colorbar(rh850contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("850 hPa Bağıl Nem (%)",labelpad=-40)



    # Set the map bounds, figurde eksenlerin sinirlari ve degerlerini alabilmek icin yazildi.
    ax.set_xlim(cartopy_xlim(rh850))
    ax.set_ylim(cartopy_ylim(rh850))
    plt.barbs(to_np(lons[::12, ::12]), to_np(lats[::12, ::12]),
                to_np(u850[0][::12, ::12]), to_np(v850[0][::12, ::12]),transform=crs.PlateCarree(),zorder=2,
                 length=5)

    plt.title("850 hPa Bağıl Nem (%)\n850 hPa Rüzgar (kt)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
    plt.savefig('rh850_wind850_{}.png'.format(time),dpi=300,transparent=False)

def rh850_t850(ncfiledir,fnamedir,timeg):
    from netCDF4 import Dataset
    from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature
    from cartopy.io.shapereader import Reader
    from matplotlib.offsetbox import  OffsetImage
    import matplotlib.image as image

    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)

    time=g_times.get_times(ncfile,timeidx=timeg,meta=False)
    time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

    baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
    baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')

    # Extract the pressure, geopotential height, and wind variables
    # Burada istediğimiz değişkenleri çekebiliriz, tablosu var.
    rh = getvar(ncfile,"rh",timeidx=timeg)
    t = getvar(ncfile,'tc',timeidx=timeg)

    rh850=vinterp(ncfile,
                field=rh,
                vert_coord="p",
                interp_levels=[850],
                extrapolate=True,
                log_p=True)

    t850=vinterp(ncfile,
                field=t,
                vert_coord="p",
                interp_levels=[850],
                extrapolate=True,
                field_type="tc",
                log_p=True)

    smooth_t850 = smooth2d(t850, 100, cenweight=4)

    # Get the lat/lon coordinatesç
    # Verinin enlem ve boylamini otomatik almasi icin yazildi.
    lats, lons = latlon_coords(t850)

    # Get the map projection information, verinin projeksiyon methodunu belirlemek icin yazildi.
    cart_proj = get_cartopy(t850)

    # Create the figure, olusturulacak goruntunun boyut ayari.
    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)
    # Download and add the states and coastlines, icindeki parametreler degisince harita cizimi degisir.
    states = NaturalEarthFeature(category="cultural", scale="10m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")
    # Yukarida belirlenen ayarlari figure eklemek icin add_feature kullanildi.
    ax.add_feature(states, linewidth=0.8, edgecolor="dimgray")
    # Sahil kenarlarini cizmeye yarar.
    ax.add_feature(COASTLINE,linewidth=0.9,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.5)
    ax.add_feature(shape_feature)
    levels = np.arange(-36,41,2)
    t850contour= plt.contour(to_np(lons), to_np(lats),to_np(smooth_t850[0]),cmap="turbo",
                transform=crs.PlateCarree(),levels=levels)
    plt.clabel(t850contour,inline=1,fontsize=8,fmt="%i",colors="black")


    rh850contour=plt.contourf(to_np(lons), to_np(lats), to_np(rh850[0]),transform=crs.PlateCarree(),cmap="Purples",levels=np.arange(0,101,10),extend="both")
    cbar= plt.colorbar(rh850contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("850 hPa Bağıl Nem (%)",labelpad=-40)



    # Set the map bounds, figurde eksenlerin sinirlari ve degerlerini alabilmek icin yazildi.
    ax.set_xlim(cartopy_xlim(t850))
    ax.set_ylim(cartopy_ylim(t850))

    plt.title("850 hPa Bağıl Nem (%) \n850 hPa Sıcaklık (°C)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
    plt.savefig('rh850_t850_{}.png'.format(time),dpi=300,transparent=False)



# for i in range(0,24):
#     wind_t2(ncfiledir,fnamedir,13)
#     vor_500h(ncfiledir,fnamedir,13)
#     slp_500h(ncfiledir,fnamedir,13)
#     slp_wind_t2(ncfiledir,fnamedir,13)
#     rh850_t850(ncfiledir,fnamedir,13)
#     rh850_wind850(ncfiledir,fnamedir,13)

rh850_t850(ncfiledir,fnamedir,5)