import pickle
import numpy as np


def read_pkl(path):
    with open(path, 'rb') as f:
        d = pickle.load(f)
    return d

def pickle_dump(obj, path):
    with open(path, mode='wb') as f:
        pickle.dump(obj, f)

def get_data(mc_data, mask_L_az, mask_L_el):
    az = mc_data["az"]
    el = mc_data["el"]
    moon_az_cond = np.abs(az-180) < mask_L_az
    moon_el_cond = np.abs(el-70) < mask_L_el
    moon_cond = moon_az_cond & moon_el_cond
    fitdata = {}
    fitdata['utime'] = mc_data['utime'][moon_cond]
    fitdata['az'] = az[moon_cond]
    fitdata['el'] = el[moon_cond]
    fitdata['phase'] = mc_data['phase'][moon_cond]
    return fitdata

def RMS(data):
  N_rms = np.sqrt(np.sum((data)**2)/data.size)
  return N_rms
