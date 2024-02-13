from lmfit import Model
import numpy as np
from tool.search_center import s_center
import pickle

def pickle_dump(obj, path):
    with open(path, mode='wb') as f:
        pickle.dump(obj, f)


class make_data(s_center):

    def pickle_dump(self, obj, path):
        with open(path, mode='wb') as f:
            pickle.dump(obj, f)

    def _get_data_(self, mc_data, mask_L_az, mask_L_el):
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


    def _get_maskd_data_(self, mc_data, mask_L_az, mask_L_el):
        az = mc_data["az"]
        el = mc_data["el"]
        moon_az_cond = np.abs(az-180) > mask_L_az
        moon_el_cond = np.abs(el-70) > mask_L_el
        moon_cond = moon_az_cond | moon_el_cond
        fitdata = {}
        fitdata['utime'] = mc_data['utime'][moon_cond]
        fitdata['az'] = az[moon_cond]
        fitdata['el'] = el[moon_cond]
        fitdata['phase'] = mc_data['phase'][moon_cond]
        return fitdata


    def _get_center_data(self, data, Laz, Lel, fitpra):
        get_data = self._get_data_(data, Laz,Lel)
        get_data["az"] = get_data["az"] - (fitpra['cen_az']-180)
        get_data["el"] = get_data["el"] - (fitpra['cen_el']-70)
        return get_data



    def _centerdata1_(self, data11):
        fitcenter = s_center()
        cen_az = data11['az'][data11['phase'].argmax()]
        cen_el = data11['el'][data11['phase'].argmax()]
        sig_az = 0.25
        sig_el = 0.25
        height = data11['phase'].max()
        fitcenter.set_hint(cen_az, cen_el, sig_az, sig_el, height)
        az_el = data11['az'], data11['el']
        pre = fitcenter.step2(data11['phase'], az_el)
        cen_az = pre['cen_az']
        cen_el = pre['cen_el']
        sig_az = pre['sig_az']
        sig_el = pre['sig_el']
        height = pre['height']
        fitcenter.set_hint(cen_az, cen_el, sig_az, sig_el, height)
        return fitcenter.step2(data11['phase'], az_el)


    def _basefunc_(self, x,c=0,d=0,a0=0,p0=0):
        ret = c + d*x
        ret += a0*np.sin(x*(2*np.pi/3)-p0)
        return ret

    def rmnoise(self, spdata, sputime, sutime):
        try:
            c_ini = spdata[0]
            d_ini = (spdata.max()-spdata.min())/(sputime.max()-sputime.min())
            result = self.model.fit(data = spdata, x = sputime, c = c_ini, d = d_ini)
            param_fit = result.best_values
            reresult = self._basefunc_(x = sutime,c=param_fit['c'],d=param_fit['d'],a0=param_fit['a0'],p0=param_fit['p0'])
        except:
            reresult = np.zeros(len(sutime))
        return reresult

    def no_noisedata(self, data, region_az, region_el, mask_az, mask_el , savename):
        # az = np.where(data["az"] < 180, data["az"]+180, data["az"]-180)
        # data["az"] = 180 + (az-180) * np.cos(np.deg2rad(data['el']))
        data11 = self._get_data_(data, 1, 1)
        cea1para = self._centerdata1_(data11)
        region_data = self._get_center_data(data, region_az, region_el, cea1para)
        hole_data = self._get_maskd_data_(region_data, mask_az, mask_el)
        diff_index_hole = np.append(-1, np.where(np.diff(hole_data["az"]) > 0))
        diff_index_hole = np.append(diff_index_hole,len(hole_data["az"])-1)
        diff_index_region = np.append(-1, np.where(np.diff(region_data["az"]) > 0))
        diff_index_region = np.append(diff_index_region,len(region_data["az"])-1)
        self.model = Model(self._basefunc_)
        self.model.make_params(verbose=True)
        self.model.set_param_hint('c')
        self.model.set_param_hint('d')
        phase_off = np.zeros(0)
        for i in range(int(len(diff_index_hole)-1)):
            phase_off = np.append(phase_off,self.rmnoise(spdata = hole_data['phase'][int(diff_index_hole[i]+1):int(diff_index_hole[i+1]+1)], sputime = hole_data['utime'][int(diff_index_hole[i]+1):int(diff_index_hole[i+1]+1)], sutime = region_data['utime'][int(diff_index_region[i]+1):int(diff_index_region[i+1]+1)]))
        dict = region_data
        dict['noise'] = phase_off
        pickle_dump(dict,savename)
