# Fonction traitant données 
# 27/03/2019 M.G.


import numpy as np
from datetime import datetime

################################################################
#			Fonction SAR			       #
################################################################


def data_sar(lon_sar,lat_sar,mask,ws_sar,angle,path_sar):
    name_sar = path_sar.split('/')[-1]
    sat_sar = name_sar.split('-')[0]
    x = np.ma.array(lon_sar)
    x.mask = np.isnan(x)    
    x[x<0] = x[x<0]+360
    y = np.ma.array(lat_sar)
    y.mask =  np.isnan(y)
    z = np.ma.masked_where(mask!=0 , ws_sar)
    w = np.ma.masked_where(mask!=0 , angle)
    xlim = (x.min(), x.max())
    ylim = (y.min(), y.max())

    day_sar = datetime.strptime(name_sar.split('-')[4][0:8],'%Y%m%d')
    day_sar_1 = day_sar.strftime('%Y_%m_%d')    
    day_sar = day_sar.strftime('%Y%m%d')
    date_sar = datetime.strptime(name_sar.split('-')[4][0:8]+name_sar.split('-')[4][9:15],'%Y%m%d%H%M%S')
    hour_sar1 = datetime.strptime(name_sar.split('-')[4][9:15],'%H%M%S')
    hour_sar_debut = hour_sar1.strftime('%H:%M:%S')
    hour_sar2 = datetime.strptime(name_sar.split('-')[5][9:15],'%H%M%S')
    hour_sar_fin = hour_sar2.strftime('%H:%M:%S')
    return name_sar,sat_sar,x,y,z,w,xlim,ylim,day_sar,day_sar_1,date_sar,hour_sar_debut,hour_sar_fin


################################################################
#			Fonction SMOS			       #
################################################################


def data_smos(x,y,smos,date_sar):
    a = np.array(smos['lon'])
    b = np.array(smos['lat'])
    e = np.where((a>=x.min()) * ( a<=x.max()) * (b>=y.min()) * (b<=y.max()))
    b1=b[e[0].min():e[0].max(),e[1].min():e[1].max()]
    a1=a[e[0].min():e[0].max(),e[1].min():e[1].max()]

    dt_mo = -999
    try:
        date_mo = np.ma.array(smos['date_mo'])
        date_mo.mask=np.isnan(date_mo)
        date_mo_smos = date_mo[e[0].min():e[0].max(),e[1].min():e[1].max()]
        moy_mo = date_mo_smos.mean() 
        date_mo_smos_h = int((moy_mo - int(moy_mo))*24)
        date_mo_smos_min = int(((moy_mo - int(moy_mo))*24 - date_mo_smos_h)*60)
        date_mo_smos_moy = datetime(smos['date_smos'].year,smos['date_smos'].month, smos['date_smos'].day, date_mo_smos_h,date_mo_smos_min)
        dt1 = date_sar-date_mo_smos_moy
        dt_mo = abs(dt1.total_seconds())
        dist_mo = smos['dist_mo'][e[0].min():e[0].max(),e[1].min():e[1].max()]
        dist_mo = dist_mo.mean()
    except:
        pass

    dt_ev = -999
    try:
        date_ev = np.ma.array(smos['date_ev'])
        date_ev.mask=np.isnan(date_ev)
        date_ev_smos = date_ev[e[0].min():e[0].max(),e[1].min():e[1].max()]
        moy_ev = date_ev_smos.mean() 
        date_ev_smos_h = int((moy_ev - int(moy_ev))*24)
        date_ev_smos_min = int(((moy_ev - int(moy_ev))*24 - date_ev_smos_h)*60)
        date_ev_smos_moy = datetime(smos['date_smos'].year,smos['date_smos'].month, smos['date_smos'].day, date_ev_smos_h,date_ev_smos_min)
        dt2 = date_sar-date_ev_smos_moy
        dt_ev = abs(dt2.total_seconds())
        dist_ev = smos['dist_ev'][e[0].min():e[0].max(),e[1].min():e[1].max()]
        dist_ev = dist_ev.mean()     
    except:
        pass

    if dt_mo != -999 and dt_ev != -999:
        if dt_mo < dt_ev:
            vent_smos = smos['spd_mo']
            dist = dist_mo
            dt = dt_mo 
            date1 = date_mo_smos_moy
        else:
            vent_smos = smos['spd_ev']
            dt = dt_ev
            dist = dist_ev
            date1 = date_ev_smos_moy
    elif dt_mo == -999 and dt_ev != -999:
        vent_smos = smos['spd_ev']
        dt = dt_ev
        dist = dist_ev
        date1 = date_ev_smos_moy
    elif dt_mo != -999 and dt_ev == -999:
        vent_smos = smos['spd_mo']
        dt = dt_mo
        dist = dist_mo 
        date1 = date_mo_smos_moy
    else:
        print('pas de colocalisation')
    vent_smos_temp = np.ma.array(vent_smos)
    vent_smos_temp.mask = np.isnan(vent_smos_temp)

    return date1,a1, b1, vent_smos_temp, dt,dist


################################################################
#			Fonction SMAP REMSS		       #
################################################################


def data_smap_remss(lon_smap_remss,lat_smap_remss,ws_smap_remss, time_smap_remss, x, y, year_smap, month_smap, day_smap, date_sar):
	
	a = np.array(lon_smap_remss)
	b = np.array(lat_smap_remss)
	lon_smap,lat_smap = np.meshgrid(a,b) 
	c = np.ma.masked_where(ws_smap_remss<0,ws_smap_remss)
	d = np.ma.masked_where(time_smap_remss<0,time_smap_remss)

	e = np.ma.where((lon_smap>=x.min()) * ( lon_smap<=x.max()) * (lat_smap>=y.min()) * (lat_smap<=y.max()))
	b1=lat_smap[e[0].min():e[0].max(),e[1].min():e[1].max()]
	a1=lon_smap[e[0].min():e[0].max(),e[1].min():e[1].max()]
	vent_smap = c[e[0].min():e[0].max(),e[1].min():e[1].max()]
	time_smap = d[e[0].min():e[0].max(),e[1].min():e[1].max()]
	dt_mo = -999
	try:
		hour_mo = np.ma.array(time_smap[:,:,1])
		hour_mo_smap = hour_mo.mean()
		h_mo_smap = int(int(hour_mo_smap)/60)
		mn_mo_smap = int((hour_mo_smap-int(hour_mo_smap))*60)
		date_mo_smap = datetime(int(year_smap),int(month_smap),int(day_smap),h_mo_smap,mn_mo_smap)
		dt1 = date_sar-date_mo_smap
		dt_mo = abs(dt1.total_seconds()) 
	except:
		pass

	dt_ev = -999
	try:
		hour_ev = np.ma.array(time_smap[:,:,0])
		hour_ev_smap = hour_ev.mean()
		h_ev_smap = int(int(hour_ev_smap)/60)
		mn_ev_smap = int((hour_ev_smap-int(hour_ev_smap))*60)
		date_ev_smap = datetime(int(year_smap),int(month_smap),int(day_smap),h_ev_smap,mn_ev_smap)
		dt2 = date_sar-date_ev_smap
		dt_ev = abs(dt2.total_seconds()) 
	except:
		pass
	if dt_mo != -999 and dt_ev != -999:
		if dt_mo < dt_ev:
			vent_smap = ws_smap_remss[:,:,1]
			dt = dt_mo 
			date1 = date_mo_smap
		else:
			vent_smap = ws_smap_remss[:,:,0]
			dt = dt_ev
			date1 = date_ev_smap
	elif dt_mo == -999 and dt_ev != -999:
		vent_smap = ws_smap_remss[:,:,0]
		dt = dt_ev
		date1 = date_ev_smap
	elif dt_mo != -999 and dt_ev == -999:
		vent_smap = ws_smap_remss[:,:,1]
		dt = dt_mo
		date1 = date_mo_smap
	else:
		print('dt_mo = -999 & dt_ev = -999')
	return(vent_smap, dt, a,b,date1, a1, b1)

