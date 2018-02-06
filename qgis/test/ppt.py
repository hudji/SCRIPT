day_of_year=230
y=((2 * math.pi) / 366) * ((day_of_year - 1) + ((10.00944988388889 - 12) / 24) ) #radians
gamma=np.radians(y) #radians

#sun declination angle
decl=0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma) - 0.006758 * np.cos (2 * gamma)\
     + 0.000907 * np.sin (2 * gamma) - 0.002697 * np.cos (3 * gamma) + 0.00148 * np.sin (2 * gamma) #radians
decl_deg= (360 / (2 * math.pi)) * decl

#eqtime
eqtime = 229.18 * (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma) - 0.014615 * np.cos(2 * gamma) - 0.040849 * np.sin(2 * gamma))  # minutes
timeoff= eqtime - 4 * long_n + 60 * 7 #minutes
tst=10 * 60 + 0 + 34.019582 / 60 + timeoff #minutes
ha=(tst /4)-180 #degree
#sun zenith angle
zenit1 =np.sin(np.deg2rad(lat_n))* np.sin(np.deg2rad(decl_deg)) + np.cos (np.deg2rad(lat_n))* np.cos(np.deg2rad(decl_deg)) * np.cos(np.deg2rad(ha))
zenit2=np.arccos(zenit1) #radians
zenit_angle= np.rad2deg(zenit2) #degrees
#sun azimuth angle
theta1= -1 * ((np.sin(np.deg2rad(lat_n)) * np.cos(np.deg2rad(zenit_angle))- np.sin(np.deg2rad(decl_deg)))/(np.cos (np.deg2rad(lat_n)) * np.sin (np.deg2rad(zenit_angle))))
theta2=np.arccos(theta1) #radians
theta3=np.rad2deg(theta2)#degree
azimuth_angle=180 - theta3 #degrees

# IC calculation
delta=azimuth_angle - aspect
IC=(np.cos(np.deg2rad(zenit_angle))* np.cos (np.deg2rad(slope))) + (np.sin(np.deg2rad(zenit_angle)) * np.sin (np.deg2rad(slope)) * np.cos(np.deg2rad(delta)))#radians
IC_D= np.rad2deg(IC)
