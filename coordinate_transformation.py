from scipy.interpolate import interp1d
import numpy as np
import math



def co_tr(data_az, data_el, data_utime, data_intensity, azoff, eloff, reference_time, reference_az, reference_el,  referencepos_el):
  az_func = interp1d(np.append(np.append(reference_time.min()-1,reference_time),reference_time.max()+1), np.append(np.append(reference_az[0],reference_az),reference_az[-1]), kind='linear')
  el_func = interp1d(np.append(np.append(reference_time.min()-1,reference_time),reference_time.max()+1),np.append(np.append(reference_el[0],reference_el),reference_el[-1]), kind='linear')
  moon_cent_az_1 = 360 - data_az - az_func(data_utime) - azoff
  moon_cent_az = np.where(moon_cent_az_1 < 0, 360 + moon_cent_az_1, moon_cent_az_1)
  EL = np.zeros(len(moon_cent_az))
  AZ = np.zeros(len(moon_cent_az))
  el = data_el - eloff
  x = np.sin(np.deg2rad(90-el)) * np.cos(np.deg2rad(-moon_cent_az))
  y = np.sin(np.deg2rad(90-el)) * np.sin(np.deg2rad(-moon_cent_az))
  z = np.cos(np.deg2rad(90-el))
  theta = np.deg2rad(el_func(data_utime)-referencepos_el)
  X = np.cos(theta)*x+np.sin(theta)*z
  Z = -np.sin(theta)*x + np.cos(theta)*z
  for i in range(len(moon_cent_az)):
      EL[i] = 90 - np.rad2deg(math.acos(Z[i]))
      AZ[i] = np.rad2deg(math.atan2(y[i], X[i]))
  AZ = np.where(AZ < 0, 360+AZ, AZ)
  return {"az": AZ, "el": EL, "eloff": eloff, "azoff": azoff,"utime": data_utime, "phase": data_intensity}
