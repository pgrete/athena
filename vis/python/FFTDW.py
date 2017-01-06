import matplotlib
matplotlib.use('Agg')

import numpy as np

import h5py
import progressbar
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob


# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

# The full quantities' names are:
#['rho' 'vel1' 'vel2' 'vel3' 'pgas' 'B1' 'B2' 'B3' 'Er' 'Fr1' 'Fr2' 'Fr3'
# 'Pr11' 'Pr22' 'Pr33' 'Pr12' 'Pr13' 'Pr23' 'Er0' 'Fr01' 'Fr02' 'Fr03'
# 'Sigma_s' 'Sigma_a']




crat=5694.76
prat=43.6907

files=sorted(glob.glob('mid*athdf'))

num_file=len(files)


nrloc=4
rloc=np.zeros((nrloc,2))
rloc[0,0]=5.0
rloc[1,0]=10.0
rloc[2,0]=15.0
rloc[3,0]=20.0

time=np.zeros(num_file)


first=0
count=0

for filename in files:
  print filename
  f=h5py.File(filename, 'r')

  x1v=f['x1v'].value
  x3v=f['x3v'].value
  
  time[count]=(f['Time'].value)*crat
  
  if first==0:
    nr=np.size(x1v)
    nphi=np.size(x3v)
    np.savetxt('phi.txt',x3v)
    np.savetxt('radius.txt',x1v)
    mode=np.zeros(nphi/2)
    for i in range(0,nphi/2):
      mode[i]=i
    np.savetxt('mode_num.txt',mode)
    for i in range(0,nrloc):
      rloc[i,1]=np.abs(x1v - rloc[i,0]).argmin()
    ST_rho_r1=np.zeros((nphi/2,num_file))
    ST_rho_r2=np.zeros((nphi/2,num_file))
    ST_rho_r3=np.zeros((nphi/2,num_file))
    ST_rho_r4=np.zeros((nphi/2,num_file))
    ST_er_r1=np.zeros((nphi/2,num_file))
    ST_er_r2=np.zeros((nphi/2,num_file))
    ST_er_r3=np.zeros((nphi/2,num_file))
    ST_er_r4=np.zeros((nphi/2,num_file))

    ST_rho_m1=np.zeros((nr,num_file))
    ST_rho_m2=np.zeros((nr,num_file))
    ST_er_m1=np.zeros((nr,num_file))
    ST_er_m2=np.zeros((nr,num_file))


    fft_er=np.zeros((nphi,nr))
    fft_er_img=np.zeros((nphi,nr))
    fft_rho=np.zeros((nphi,nr))
    fft_rho_img=np.zeros((nphi,nr))
    first=1
#calculate the momentum gradient
  rho= f['rho'].value
  Er=f['Er'].value
#calculate FFT for each radius

  for i in range(nr):
    fft_temp=np.fft.fft(rho[:,i]-np.average(rho[:,i]))
    fft_rho[:,i]=np.real(fft_temp)
    fft_rho_img[:,i]=np.imag(fft_temp)
    fft_temp=np.fft.fft(Er[:,i]-np.average(Er[:,i]))
    fft_er[:,i]=np.real(fft_temp)
    fft_er_img[:,i]=np.imag(fft_temp)

  absfft_rho=(fft_rho**2+fft_rho_img**2)**0.5
  absfft_er=(fft_er**2+fft_er_img**2)**0.5

  ST_rho_r1[:,count]=absfft_rho[0:nphi/2,rloc[0,1]]
  ST_rho_r2[:,count]=absfft_rho[0:nphi/2,rloc[1,1]]
  ST_rho_r3[:,count]=absfft_rho[0:nphi/2,rloc[2,1]]
  ST_rho_r4[:,count]=absfft_rho[0:nphi/2,rloc[3,1]]

  ST_er_r1[:,count]=absfft_er[0:nphi/2,rloc[0,1]]
  ST_er_r2[:,count]=absfft_er[0:nphi/2,rloc[1,1]]
  ST_er_r3[:,count]=absfft_er[0:nphi/2,rloc[2,1]]
  ST_er_r4[:,count]=absfft_er[0:nphi/2,rloc[3,1]]

  ST_rho_m1[:,count]=absfft_rho[1,:]
  ST_rho_m2[:,count]=absfft_rho[2,:]

  ST_er_m1[:,count]=absfft_rho[1,:]
  ST_er_m2[:,count]=absfft_rho[2,:]
  count=count+1

np.savetxt('FFT_rho_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_rho_r1,fmt='%4.5e')
np.savetxt('FFT_rho_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_rho_r2,fmt='%4.5e')
np.savetxt('FFT_rho_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_rho_r3,fmt='%4.5e')
np.savetxt('FFT_rho_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_rho_r4,fmt='%4.5e')

np.savetxt('FFT_er_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_er_r1,fmt='%4.5e')
np.savetxt('FFT_er_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_er_r2,fmt='%4.5e')
np.savetxt('FFT_er_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_er_r3,fmt='%4.5e')
np.savetxt('FFT_er_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_er_r4,fmt='%4.5e')

np.savetxt('FFT_rho_m1.txt',ST_rho_m1,fmt='%4.5e')
np.savetxt('FFT_rho_m2.txt',ST_rho_m2,fmt='%4.5e')

np.savetxt('time.txt',time)
