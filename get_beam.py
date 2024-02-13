from tool.noise_remove import make_data
from scipy import interpolate
import numpy as np

class beam_data(make_data):
  def info(self, N,Laz,Lel):
    self.make_data = make_data()
    self.N = N/2
    self.Laz = Laz
    self.Lel = Lel
    self.naz = np.linspace(self.N-1, -self.N, int(2*self.N))
    self.nel = np.linspace(self.N-1, -self.N, int(2*self.N))
    self.delaz = self.Laz/self.N
    self.delel = self.Lel/self.N
    self.az = 180 + self.naz*self.delaz
    self.el = 70 + self.nel*self.delel
    self.az_grid, self.el_grid = np.meshgrid(self.az,self.el)


  def inp(self, nrm_data, saverbf, kernel):
    nrm_data['phase'] = nrm_data['phase'] - nrm_data['noise']
    data11 = self._get_data_(nrm_data, 1, 1)
    cea1para = self._centerdata1_(data11)
    self.data = self._get_center_data(nrm_data, self.Laz+0.1, self.Lel+0.1, cea1para)
    rbf_l = interpolate.Rbf(self.data['az'], self.data['el'], self.data['phase'], function = kernel)
    zz_new = rbf_l(self.az_grid, self.el_grid)
    dict = {}
    dict["az_grid"] = self.az_grid
    dict["el_grid"] = self.el_grid
    dict["Rbf"] = zz_new
    self.pickle_dump(dict, saverbf)
    return zz_new


  def inp2(self, nrm_data, saverbf, kernel):
    nrm_data['phase'] = nrm_data['phase'] - nrm_data['noise']
    data11 = self._get_data_(nrm_data, 1, 1)
    cea1para = self._centerdata1_(data11)
    self.data = self._get_center_data(nrm_data, self.Laz+0.1, self.Lel+0.1, cea1para)
    rbf_l = interpolate.Rbf(self.data['az'], self.data['el'], self.data['phase'], function= kernel)
    zz_new = rbf_l(self.az_grid, self.el_grid)
    dict = {}
    dict["az_grid"] = self.az_grid
    dict["el_grid"] = self.el_grid
    dict["Rbf"] = zz_new
    self.pickle_dump(dict, saverbf)
    return zz_new


  def _dft_(self, saz, sel, fase_arr):
    e = np.outer(np.exp(-1j*(self.nel*np.pi*sel)/self.N), np.exp(-1j*(self.naz*np.pi*saz)/self.N))
    f = np.sum(fase_arr * e)
    return f*(self.Laz*self.Lel)/(self.N * self.N)

  def get_DFT(self, phase):
    dft_arr = np.zeros((int(2*self.N), int(2*self.N)), dtype=complex)
    for i in range(int(2*self.N)):
        for j in range(int(2*self.N)):
            dft_arr[i, j] = self._dft_(self.N-1-j, self.N-1-i, phase)
    return dft_arr

  def _idft_(self, saz, sel, kernel, kil):
    naz_grid, nel_grid = np.meshgrid(self.naz, self.nel)
    e = np.outer(np.exp(1j*((self.nel*np.pi*sel)/self.N)),np.exp(1j*((self.naz*np.pi*saz)/self.N)))
    f = np.sum(kernel[(np.abs(nel_grid) < kil)&(np.abs(naz_grid) < kil)]*e[(np.abs(nel_grid) < kil)&(np.abs(naz_grid) < kil)])
    return f/(4*self.Laz*self.Lel)

  def get_IDFT(self, kernel, kil):
    K_r_a = np.zeros((int(2*self.N), int(2*self.N)), dtype=complex)
    for i in range(int(2*self.N)):
      for j in range(int(2*self.N)):
        K_r_a[i, j]=self._idft_(-self.N+1+j, -self.N+1+i, kernel, kil)
    return K_r_a

  def Epara(self, theta):
    return 2*1.8*np.sqrt(1.8**2-np.sin(theta)**2)/(np.sqrt(1.8**2-np.sin(theta)**2)+1.8**2*np.cos(theta))

  def Ever(self, theta):
    return 2*np.sqrt(1.8**2-np.sin(theta)**2)/(np.sqrt(1.8**2-np.sin(theta)**2)+np.cos(theta))

  def I(self, r):
    theta = np.arcsin(r)
    return ((self.Epara(theta)**2+self.Ever(theta)**2)/2)*np.cos(theta)/np.sqrt(1.8**2-np.sin(theta)**2)

  def moon(self):
    r = np.sqrt((self.az_grid-180)**2+(self.el_grid-70)**2)
    r_R = np.where(r > 0.25, -1, r/0.25)
    MOON = self.I(r_R)
    return self.get_DFT(MOON)
