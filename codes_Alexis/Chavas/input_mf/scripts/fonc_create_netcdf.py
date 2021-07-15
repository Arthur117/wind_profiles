# Script permettant la création d'un fichier netcdf intermédiare contenant les données nécessaires pour 
# faire une étude statistique

def create_netcdf_smos(file_net,name_cycl, name_sar, name_smos, hour_sar_debut, hour_sar_fin, date_smos_moyenne, lon_smosss, lat_smosss, vent_sarmoy, vent_sarmed, vent_sarstd, ws_smosss, n_sarvalid, dt, angle_sarmoy,sat,cat,sat_sar,sat_rad ):
	file_net.description = 'Data processed from Sentinel-1 to SMOS resolution for tropical hurricane '+name_cycl
	file_net.source_sar = name_sar
	file_net.source_smos = name_smos
	file_net.cyclone = name_cycl
	file_net.hour_sar = 'Begin : '+hour_sar_debut+' End : '+hour_sar_fin
	file_net.hour_smos = date_smos_moyenne
	file_net.category = str(cat)
	file_net.SAR = sat+' '+sat_sar
	file_net.SMAP = sat_rad

	file_net.createDimension('dim1', len(lon_smosss[0,:]))
	file_net.createDimension('dim2', len(lat_smosss[:,0]))
	file_net.createDimension('date', None)
	

	longitude = file_net.createVariable('Lon', 'f4', 'dim1')
	longitude.units = "degrees_east"
	longitude[:] = lon_smosss[0,:]

	latitude = file_net.createVariable('Lat', 'f4', 'dim2')  
	latitude.units = "degrees_north"
	latitude[:] = lat_smosss[:,0]

	ventsarmoy = file_net.createVariable('WindSpeed_sar_mean', 'f4', ('dim2','dim1'))
	ventsarmoy.units = "m/s"
	ventsarmoy.long_name = " Mean value of wind speed "
	ventsarmoy[:]=vent_sarmoy

	mediane = file_net.createVariable('WindSpeed_sar_median', 'f4', ('dim2','dim1') )
	mediane.units = "m/s"
	mediane.long_name = " Median value of wind speed "
	mediane[:]=vent_sarmed

	ecart_type = file_net.createVariable('WindSpeed_sar_std', 'f4', ('dim2','dim1'))
	ecart_type.units = "m/s"
	ecart_type.long_name = " Standard deviation "
	ecart_type[:] = vent_sarstd

	ventsmos = file_net.createVariable('WindSpeed_smos', 'f4', ('dim2','dim1'))
	ventsmos.units = " m/s "
	ventsmos.long_name = " Wind speed by smos acquisiion "
	ventsmos[:] = ws_smosss

	nsarvalid = file_net.createVariable('number_points', 'f4',  ('dim2','dim1'))
	nsarvalid.units = " no dimension "
	nsarvalid.long_name = " Number of points "
	nsarvalid[:] = n_sarvalid

	date = file_net.createVariable('delta_time','f4','date')
	date.units = "s"
	date.long_name = "Time difference between sar and smos acquisition"
	date[:] = dt

	anglesarmoy = file_net.createVariable('sar_mean_incidence_angle', 'f4',  ('dim2','dim1'))
	anglesarmoy.units = "degrees"
	anglesarmoy.long_name = " Incidence angle "
	anglesarmoy[:]=angle_sarmoy

	return

def create_netcdf_smap(file_net,name_cycl, name_sar, name_smap, hour_sar_debut, hour_sar_fin, date_smap_moyenne, lon_smapss, lat_smapss, vent_sarmoy, vent_sarmed, vent_sarstd, ws_smapss, n_sarvalid, dt, angle_sarmoy,sat,cat,sat_sar,sat_rad ):
	file_net.description = 'Data processed from Sentinel-1 to SMAP resolution for tropical hurricane '+name_cycl
	file_net.source_sar = name_sar
	file_net.source_smos = name_smap
	file_net.cyclone = name_cycl
	file_net.hour_sar = 'Begin : '+hour_sar_debut+' End : '+hour_sar_fin
	file_net.hour_smap = date_smap_moyenne
	file_net.category = str(cat)
	file_net.SAR = sat+' '+sat_sar
	file_net.SMAP = sat_rad

	file_net.createDimension('dim1', len(lon_smapss[0,:]))
	file_net.createDimension('dim2', len(lat_smapss[:,0]))
	file_net.createDimension('date', None)

	longitude = file_net.createVariable('Lon', 'f4', 'dim1')
	longitude.units = "degrees_east"
	longitude[:] = lon_smapss[0,:]

	latitude = file_net.createVariable('Lat', 'f4', 'dim2')  
	latitude.units = "degrees_north"
	latitude[:] = lat_smapss[:,0]

	ventsarmoy = file_net.createVariable('WindSpeed_sar_mean', 'f4', ('dim2','dim1'))
	ventsarmoy.units = "m/s"
	ventsarmoy.long_name = " Mean value of wind speed "
	ventsarmoy[:]=vent_sarmoy

	mediane = file_net.createVariable('WindSpeed_sar_median', 'f4',('dim2','dim1') )
	mediane.units = "m/s"
	mediane.long_name = " Median value of wind speed "
	mediane[:]=vent_sarmed

	ecart_type = file_net.createVariable('WindSpeed_sar_std', 'f4', ('dim2','dim1'))
	ecart_type.units = "m/s"
	ecart_type.long_name = " Standard deviation "
	ecart_type[:] = vent_sarstd

	ventsmap = file_net.createVariable('WindSpeed_smap', 'f4', ('dim2','dim1'))
	ventsmap.units = " m/s "
	ventsmap.long_name = " Wind speed by smap acquisiion "
	ventsmap[:] = ws_smapss

	nsarvalid = file_net.createVariable('number_points', 'f4', ('dim2','dim1'))
	nsarvalid.units = " no dimension "
	nsarvalid.long_name = " Number of points "
	nsarvalid[:] = n_sarvalid

	date = file_net.createVariable('delta_time','f4','date')
	date.units = "s"
	date.long_name = "Time difference between sar and smap acquisition"
	date[:] = dt

	anglesarmoy = file_net.createVariable('sar_mean_incidence_angle', 'f4', ('dim2','dim1'))
	anglesarmoy.units = "degrees"
	anglesarmoy.long_name = " Incidence angle "
	anglesarmoy[:]=angle_sarmoy

	return


