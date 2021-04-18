from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature,COASTLINE,ShapelyFeature,COLORS,LAND,OCEAN,BORDERS
from cartopy.io.shapereader import Reader
from matplotlib.offsetbox import  OffsetImage
import matplotlib.image as image
import datetime
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim,smooth2d,vinterp,g_times,ALL_TIMES)
# ncfiledir=input('Ncfile addres')
# fnamedir=input('fname')

ncfiledir = "/home/durmaz/wrfout/wrfout_d01_2021-02-14.nc"
fnamedir = "/home/durmaz/Visulazation/wrf_out_process/shapefiles/gadm36_TUR_1.shp"

def vor_500h(ncfiledir,fnamedir,timeg):

    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)

    p = getvar(ncfile, "pressure",timeidx=timeg)
    slp = getvar(ncfile, "slp",units='hPa',timeidx=timeg)
    z=getvar(ncfile,"z",units="m",timeidx=timeg)
    vor=getvar(ncfile,"avo",timeidx=timeg)
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
    vor500=vinterp(ncfile,
                field=vor,
                vert_coord="p",
                interp_levels=[500],
                extrapolate=True,
                log_p=True)

    smooth_ht500 = smooth2d(ht500, 2, cenweight=4)
    lats, lons = latlon_coords(slp)

    cart_proj = get_cartopy(slp)

    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("/home/durmaz/Visulazation/wrf_out_process/logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")
    ax.add_feature(states, linewidth=0.4, edgecolor="black")
    ax.add_feature(COASTLINE,linewidth=0.4,edgecolor="black")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="black",linewidth=0.4)
    ax.add_feature(shape_feature)

    contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_ht500[0]/10),
                            colors="dimgray",
                            transform=crs.PlateCarree(),linewidths=1,levels=np.arange(476,604,4))
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i",colors="black")
    levels = np.arange(-12,13,2)
    vor500contour= plt.contourf(to_np(lons), to_np(lats),to_np(vor500[0]),extend="both",cmap="coolwarm",
                transform=crs.PlateCarree(),levels=levels)
    cbar= plt.colorbar(vor500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)

    cbar.set_label("500 hPa Vortisiti (10-5 s-1)",labelpad=-40)


    ax.set_xlim(cartopy_xlim(slp))
    ax.set_ylim(cartopy_ylim(slp))

    plt.title("500 hPa Vortisiti (10-5 s-1) \n500 hPa Jeopotansiyel Yüksekliği (dam)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
    plt.savefig('vor_500h_{}.png'.format(time),dpi=300,transparent=False)

def pre_slp_wind(ncfiledir,fnamedir):
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    alltimes= getvar(ncfile, "times",ALL_TIMES)

    for timeindex in range(3,len(alltimes)-18,3):

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

        levels = [0.25, 0.50, 1, 2, 5, 10, 15, 20,25, 30, 35, 40, 45, 50]
        precip= plt.contourf(to_np(lons), to_np(lats),np.ma.masked_less((pre_total),0.2),extend="both",cmap="summer",
                    transform=crs.PlateCarree(),levels=levels,zorder=2)
        cbar= plt.colorbar(precip,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
        cbar.set_label("Yağış (mm) ",labelpad=-40)

        plt.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
        to_np(u[::10, ::10]), to_np(v[::10, ::10]),
        transform=crs.PlateCarree(), length=4,zorder=2,linewidth=0.4)

        slpcontours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                    colors="maroon", transform=crs.PlateCarree(),linewidths=0.6,zorder=2)
        plt.clabel(slpcontours, inline=1, fontsize=5, fmt="%i",colors="maroon")

        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))

        plt.title("Yağış (mm)\n10m Rüzgar\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
        plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

        plt.savefig('pre_slp_wind_{}.png'.format(time),dpi=300,transparent=False)

def pre(ncfiledir,fnamedir):
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    alltimes= getvar(ncfile, "times",ALL_TIMES)

    for timeindex in prange(3,len(alltimes),3):

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

        levels = [0.25, 0.50, 1, 2, 5, 10, 15, 20,25, 30, 35, 40, 45, 50]
        precip= plt.contourf(to_np(lons), to_np(lats),np.ma.masked_less((pre_total),0.2),extend="both",cmap="summer",
                    transform=crs.PlateCarree(),levels=levels,zorder=2)
        cbar= plt.colorbar(precip,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
        cbar.set_label("Yağış (mm) ",labelpad=-40)


        ax.set_xlim(cartopy_xlim(slp))
        ax.set_ylim(cartopy_ylim(slp))
        plt.title("Yağış (6h)(mm)",loc="left",pad=10,fontsize=12)
        plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

        plt.savefig('pre_{}.png'.format(time.replace('-','_').replace(':','_')),dpi=300,transparent=False)

def slp_500h_500t(ncfiledir,fnamedir,timeg):
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)

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

    smooth_slp = smooth2d(slp,100, cenweight=4)
    smooth_t500 = smooth2d(t500, 3, cenweight=4)

    lats, lons = latlon_coords(slp)

    cart_proj = get_cartopy(slp)

    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")
    ax.add_feature(states, linewidth=0.4, edgecolor="dimgray")
    ax.add_feature(COASTLINE,linewidth=0.5,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.2)
    ax.add_feature(shape_feature)

    contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                            colors="black",levels=levels
                            transform=crs.PlateCarree(),linewidths=1)
    plt.clabel(contours, inline=1, fontsize=10, fmt="%i",colors="black")

    levels = np.arange(476,604,4)
    ht500contour= plt.contourf(to_np(lons), to_np(lats),to_np(ht500[0]/10),extend="both",cmap="rainbow",
                transform=crs.PlateCarree(),levels=levels)
    cbar= plt.colorbar(ht500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("500 hPa Jeopotansiyel Yükseklik (metre)",labelpad=-40)

    tcontour=plt.contour(to_np(lons), to_np(lats), to_np(smooth_t500[0]),colors="black",
                            transform=crs.PlateCarree(),linewidths=0.5,linestyles="--")
    plt.clabel(tcontour, inline=1, fontsize=8, fmt="%i",colors="black")

    ax.set_xlim(cartopy_xlim(slp))
    ax.set_ylim(cartopy_ylim(slp))

    plt.title("500 hPa Sıcaklığı (°C)\n500 hPa Jeopotansiyel Yüksekliği (m)\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
    plt.savefig('slp_500h_500t_{}.png'.format(time),dpi=300,transparent=False)

def slp_500h(ncfiledir,fnamedir,timeg):
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)

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

    smooth_slp = smooth2d(slp,100, cenweight=4)
    smooth_t500 = smooth2d(t500, 3, cenweight=4)

    lats, lons = latlon_coords(slp)

    cart_proj = get_cartopy(slp)

    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("/home/durmaz/Visulazation/wrf_out_process/logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")
    ax.add_feature(states, linewidth=0.4, edgecolor="dimgray")
    ax.add_feature(COASTLINE,linewidth=0.5,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.2)
    ax.add_feature(shape_feature)

    levels = np.arange(950,6040,40)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                            colors="black",
                            transform=crs.PlateCarree(),linewidths=1,levels=np.arange(950,1100,2))
    plt.clabel(contours, inline=1, fontsize=7, fmt=lambda p: f'{str(int(p))[-2:]}')
    

    levels = np.arange(476,604,4)
    ht500contour= plt.contourf(to_np(lons), to_np(lats),to_np(ht500[0]/10),extend="both",cmap="rainbow",
                transform=crs.PlateCarree(),levels=levels)
    cbar= plt.colorbar(ht500contour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)
    cbar.set_label("500 hPa Jeopotansiyel Yükseklik (dam)",labelpad=-40)

    ax.set_xlim(cartopy_xlim(slp))
    ax.set_ylim(cartopy_ylim(slp))

    plt.title("500 hPa Jeopotansiyel Yüksekliği (m)\nDeniz Seviyesi Basıncı (mb)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)
    plt.savefig('slp_500h_{}.png'.format(time),dpi=300,transparent=False)

def slp_wind_t2(ncfiledir,fnamedir,timeg):
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    slp = getvar(ncfile, "slp",units='hPa',timeidx=timeg)
    t= getvar(ncfile,'T2',timeidx=timeg)
    t=t-273.15
    u,v= getvar(ncfile,'uvmet10',timeidx=timeg)

    time=g_times.get_times(ncfile,timeidx=timeg,meta=False)
    time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

    baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
    baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')

    smooth_slp = smooth2d(slp, 100, cenweight=4)


    lats, lons = latlon_coords(slp)
    cart_proj = get_cartopy(slp)

    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)

    states = NaturalEarthFeature(category="cultural", scale="50m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")

    ax.add_feature(states, linewidth=0.4, edgecolor="dimgray")

    ax.add_feature(COASTLINE,linewidth=0.8,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.2)
    ax.add_feature(shape_feature)

    levels = np.arange(950,1100,2)
    slpcontours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp),
                            colors="maroon",
                            transform=crs.PlateCarree(),linewidths=1,levels=levels)

    plt.clabel(slpcontours, inline=1, fontsize=10, fmt=lambda p: f'{str(int(p))[-2:]}',colors="black")
    levels = np.arange(-30,44,2)
    t2mcontour= plt.contourf(to_np(lons), to_np(lats),to_np(t),extend="both",cmap="jet",
                transform=crs.PlateCarree(),levels=levels)
    cbar= plt.colorbar(t2mcontour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)

    cbar.set_label("2 metre Sıcaklığı (°C) ",labelpad=-40)


    plt.barbs(to_np(lons[::12,::12]), to_np(lats[::12,::12]),
            to_np(u[::12, ::12]), to_np(v[::12, ::12]),
            transform=crs.PlateCarree(), length=5)

    ax.set_xlim(cartopy_xlim(slp))
    ax.set_ylim(cartopy_ylim(slp))
    plt.title("2 metre Sıcaklığı (°C)\n10 metre rüzgar (kt)\nDeniz Seviyesi Basıncı(hPa)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

    plt.savefig('slp_wind_t2_{}.png'.format(time),dpi=300,transparent=False)

def wind_t2(ncfiledir,fnamedir,timeg):
    
    ncfile = Dataset(str(ncfiledir))
    fname = str(fnamedir)
    slp = getvar(ncfile, "slp",units='hPa',timeidx=timeg)
    t= getvar(ncfile,'T2',timeidx=timeg)
    t=t-273.15
    u,v= getvar(ncfile,'uvmet10',timeidx=timeg)

    time=g_times.get_times(ncfile,timeidx=timeg,meta=False)
    time = datetime.datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%f000')

    baslangictime=g_times.get_times(ncfile,timeidx=0,meta=False)
    baslangictime = datetime.datetime.strptime(str(baslangictime),'%Y-%m-%dT%H:%M:%S.%f000')



    lats, lons = latlon_coords(slp)
    cart_proj = get_cartopy(slp)

    fig = plt.figure(figsize=(10,6))

    ax = plt.axes(projection=cart_proj)
    logo=plt.imread("logo.png")

    ax.figure.figimage(logo, 380, 400,alpha=0.7)

    states = NaturalEarthFeature(category="cultural", scale="50m",
                                facecolor="none",
                                name="admin_0_boundary_lines_land")

    ax.add_feature(states, linewidth=0.4, edgecolor="dimgray")

    ax.add_feature(COASTLINE,linewidth=0.8,edgecolor="dimgray")
    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                    crs.PlateCarree(), facecolor='none',edgecolor="dimgray",linewidth=0.2)
    ax.add_feature(shape_feature)

    levels = np.arange(-30,44,2)
    t2mcontour= plt.contourf(to_np(lons), to_np(lats),to_np(t),extend="both",cmap="jet",
                transform=crs.PlateCarree(),levels=levels)
    cbar= plt.colorbar(t2mcontour,fraction=0.071,orientation="horizontal",pad=.07,aspect=50)

    cbar.set_label("2 metre Sıcaklığı (°C) ",labelpad=-40)


    plt.barbs(to_np(lons[::12,::12]), to_np(lats[::12,::12]),
            to_np(u[::12, ::12]), to_np(v[::12, ::12]),
            transform=crs.PlateCarree(), length=5)

    ax.set_xlim(cartopy_xlim(slp))
    ax.set_ylim(cartopy_ylim(slp))
    plt.title("2 metre Sıcaklığı (°C)\n10 metre rüzgar (kt)",loc="left",pad=10,fontsize=12)
    plt.title("Başlangıç: {} UTC\nGeçerli: {} UTC".format(baslangictime,time),loc="right",pad=10,fontsize=10)

    plt.savefig('wind_t2_{}.png'.format(time),dpi=300,transparent=False)



wind_t2(ncfiledir,fnamedir,13)
vor_500h(ncfiledir,fnamedir,13)
slp_500h_500t(ncfiledir,fnamedir,13)
slp_500h(ncfiledir,fnamedir,13)
slp_wind_t2(ncfiledir,fnamedir,13)
