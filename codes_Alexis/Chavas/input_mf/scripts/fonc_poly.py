# Fonction definissant les polygones des zones cycloniques et les traçant
# 17/05/2019 M.G


from shapely.geometry.polygon import LinearRing, Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def poly():
	area_1 = Polygon([(-110,40),(-110,65),(-43,65),(-43,5),(-75,5),(-90,15),(-100,20)])
	area_2 = Polygon([(-140,65),(-110,65),(-110,40),(-100,20),(-90,15),(-75,5),(-80,5),(-140,5)])
	area_3 = Polygon([(-180, 0), (-140, 0), (-140, 65),(-180, 65)])
	area_4 = Polygon([(100, 0), (100, 60), (180, 60),(180, 0)])
	area_5 = Polygon([(50,3), (100,3),(100,60), (50,60)])
	area_6 = Polygon([(20,-5), (90,-5),(90,-30), (20,-30)])
	area_7 = Polygon([(90,0), (160,0),(160,-37), (90,-37)])
	area_11 = Polygon([(-180,0), (-120,0),(-120,-25), (-180,-25)])
	area_11b = Polygon([(180,0), (160,0),(160,-25), (180,-25)])
	area_12 = Polygon([(-180,-25), (-120,-25),(-120,-40), (-180,-40)])
	area_12b = Polygon([(180,-25), (160,-25),(160,-40), (180,-40)])

	figure(figsize = (20,12))
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.stock_img()
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
	gl.xlabels_top = False
	gl.ylabels_right = False
	gl.xformatter = LONGITUDE_FORMATTER
	gl.yformatter = LATITUDE_FORMATTER

	x_1,y_1 = area_1.exterior.xy
	plt.plot(x_1, y_1, color='k', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 1')
	x_2,y_2 = area_2.exterior.xy
	plt.plot(x_2, y_2, color='grey', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 2')
	x_3,y_3 = area_3.exterior.xy
	plt.plot(x_3, y_3, color='maroon', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 3')
	x_4,y_4 = area_4.exterior.xy
	plt.plot(x_4, y_4, color='coral', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 4')
	x_5,y_5 = area_5.exterior.xy
	plt.plot(x_5, y_5, color='g', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 5')
	x_6,y_6 = area_6.exterior.xy
	plt.plot(x_6, y_6, color='dodgerblue', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 6')
	x_7,y_7 = area_7.exterior.xy
	plt.plot(x_7, y_7, color='darkviolet', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 7,8,9,10')
	x_11,y_11 = area_11.exterior.xy
	plt.plot(x_11, y_11, color='crimson', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 11')
	x_12,y_12 = area_12.exterior.xy
	plt.plot(x_12,y_12, color='darkorange', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2,label='area 12')
	x_11b,y_11b = area_11b.exterior.xy
	plt.plot(x_11b, y_11b, color='crimson', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2)
	x_12b,y_12b = area_12b.exterior.xy
	plt.plot(x_12b,y_12b, color='darkorange', alpha=0.7,linewidth=2, solid_capstyle='round', zorder=2)
	plt.legend()
