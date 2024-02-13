"""search center by fitting 2Dgaussian"""

import numpy as np
from lmfit import Model


class s_center():
    def __init__(self):
        self.cen_az= 180
        self.cen_el=70
        self.sig_az=0.4
        self.sig_el=0.4
        self.height=0.15

    def _fitfunc_(azel, cen_az, cen_el, sig_az, sig_el, height):
        az,el = azel
        g = height*np.exp(-(((cen_az-az)/sig_az)**2 + ((cen_el-el)/sig_el)**2)/2.0)
        return g

    # def step1(self):
    #     self.model = Model(_fitfunc_)
    #     self.model.make_params(verbose=True)
    #     self.model.set_param_hint('cen_az')
    #     self.model.set_param_hint('cen_el')
    #     self.model.set_param_hint('height')
    #     self.model.set_param_hint('sig_az')
    #     self.model.set_param_hint('sig_el')
    #     self.model.set_param_hint('height')

    def set_hint(self, cen_az, cen_el, sig_az, sig_el, height):
        """input hint"""
        self.cen_az=cen_az
        self.cen_el=cen_el
        self.sig_az=sig_az
        self.sig_el=sig_el
        self.height=height


    def step2(self, intencity, az_el):
        model = Model(s_center._fitfunc_)
        model.make_params(verbose=True)
        model.set_param_hint('cen_az')
        model.set_param_hint('cen_el')
        model.set_param_hint('height')
        model.set_param_hint('sig_az')
        model.set_param_hint('sig_el')
        result = model.fit(data = intencity,  azel = az_el, cen_az = self.cen_az, cen_el = self.cen_el, sig_az = self.sig_az, sig_el = self.sig_el, height = self.height)
        return result.best_values
