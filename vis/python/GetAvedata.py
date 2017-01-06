import matplotlib
matplotlib.use('Agg')

import yt
import sys
sys.settrace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os
import glob
import h5py



# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

Crat=5694.76
Prat=43.6907
kappaes=5.04655e5
#rs=2GM/c^2==1
Medd = 20*np.pi*Crat/kappaes
Ledd = 2* np.pi * Crat * Crat * Crat/kappaes

vol_func = lambda rm,rp,thetam,thetap: \
           (1.0/3.0)*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

rarea_func = lambda r,thetam,thetap: \
           r**2.0 * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

ni=11580
no=11589


first=0

avedata={}


for fi in range(ni,no+1):
  filename='average_'+'{:05d}'.format(fi)+'.athdf'
  print filename
  f=h5py.File(filename, 'r')

  quantities=f.keys()
  
  rhovr=f['RhoVr'].value
  rhovphi=f['RhoVphi'].value
  rhovtheta=f['RhoVtheta'].value
  rho=f['rho'].value
  Ekin1=f['Ekin1'].value
  Ekin3=f['Ekin3'].value
  PB=f['PB'].value
  rhosq=f['rhosq'].value
  PBsq=f['PBsq'].value
  
  

  meanvr=rhovr/rho
  meanvphi=rhovphi/rho
  dEkin1=Ekin1-0.5*rho*meanvr*meanvr
  dEkin3=Ekin3-0.5*rho*meanvphi*meanvphi
  meanrhovphivr=meanvphi*rhovr
  meanrhovphivtheta=meanvphi*rhovtheta
  
  sigmarho=(rhosq-rho*rho)**0.5
  sigmaPB=(PBsq-PB*PB)**0.5
  meanrhosq=rho*rho
  meanpbsq=PB*PB
  meanrhoPB=rho*PB

  if first==0:
    for q in quantities:
      avedata[q]=f[q].value
    avedata['dEkin1']=dEkin1
    avedata['dEkin3']=dEkin3
    avedata['meanrhovrvphi']=meanrhovphivr
    avedata['meanrhovphivtheta']=meanrhovphivtheta
    avedata['meanvr']=meanvr
    avedata['meanvphi']=meanvphi
    avedata['meanrhosq']=meanrhosq
    avedata['meanpbsq']=meanpbsq
    avedata['meanrhoPB']=meanrhoPB
    avedata['sigmarho']=sigmarho
    avedata['sigmaPB']=sigmaPB
    first=1
  else:
    for q in quantities:
      if q != 'x1f' and q != 'x2f' and q != 'x1v' and q != 'x2v' and q != 'Time':
        avedata[q] = avedata[q] + f[q].value
    avedata['dEkin1']=avedata['dEkin1']+dEkin1
    avedata['dEkin3']=avedata['dEkin3']+dEkin3
    avedata['meanrhovrvphi']=avedata['meanrhovrvphi']+meanrhovphivr
    avedata['meanrhovphivtheta']=avedata['meanrhovphivtheta']+meanrhovphivtheta
    avedata['meanvr']=avedata['meanvr']+meanvr
    avedata['meanvphi']=avedata['meanvphi']+meanvphi
    avedata['meanrhosq']=avedata['meanrhosq']+meanrhosq
    avedata['meanpbsq']=avedata['meanpbsq']+meanpbsq
    avedata['meanrhoPB']=avedata['meanrhoPB']+meanrhoPB
    avedata['sigmarho']=avedata['sigmarho']+sigmarho
    avedata['sigmaPB']=avedata['sigmaPB']+sigmaPB

  # add items in quantities

# add to quantities
quantities.append('dEkin1')
quantities.append('dEkin3')
quantities.append('meanrhovrvphi')
quantities.append('meanrhovphivtheta')
quantities.append('meanvr')
quantities.append('meanvphi')
quantities.append('meanrhosq')
quantities.append('meanpbsq')
quantities.append('meanrhoPB')
quantities.append('sigmarho')
quantities.append('sigmaPB')

ntot=no-ni+1

for q in quantities:
   if q != 'x1f' and q != 'x2f' and q != 'x1v' and q != 'x2v' and q != 'Time':
      avedata[q] = avedata[q]/ntot


# save the time averaged data
outfile='timeavedata_'+'{:05d}'.format(ni)+'_'+'{:05d}'.format(no)+'.athdf'
outf=h5py.File(outfile, 'w')

for q in quantities:
    outf.create_dataset(q,data=avedata[q],dtype='>f8',shape=avedata[q].shape)

outf.close()


## create radial profiles

x1f=avedata['x1f']
x2f=avedata['x2f']
x1v=avedata['x1v']
x2v=avedata['x2v']

rho=avedata['rho']
PB=avedata['PB']
B1=avedata['B1']
B2=avedata['B2']
B3=avedata['B3']
Maxwell=avedata['Maxwell']
Reynolds=avedata['Reynolds']
Er=avedata['Er']
pgas=avedata['pgas']
vel3=avedata['vel3']
vel1=avedata['vel1']
Fr1=avedata['Fr1']
radacc=avedata['Radacc']*Prat
kappa_a=avedata['kappa_a']
kappa_s=avedata['kappa_s']
kappa=(kappa_a+kappa_s)/kappaes
rhov=avedata['RhoVr']
rhovin=avedata['RhoVin']
rhovout=avedata['RhoVout']
Ekin=avedata['Ekin1']+avedata['Ekin2']+avedata['Ekin3']

nr=len(x1v)
ntheta=len(x2v)

vol=np.zeros((ntheta,nr))
x1area=np.zeros((ntheta,nr+1))

for j in range(0,ntheta):
  for i in range(0,nr):
     cellsize=vol_func(x1f[i],x1f[i+1],x2f[j],x2f[j+1])
     vol[j,i]=cellsize

for j in range(0,ntheta):
  for i in range(0,nr+1):
     x1area[j,i]=rarea_func(x1f[i],x2f[j],x2f[j+1])


radialprofile=np.zeros((21,nr))
radialprofile2=np.zeros((21,nr))
radialprofile[0,:]=x1v
radialprofile2[0,:]=x1v


for i in range(0,nr):
  radialprofile[1,i]=np.average(rho[:,i],weights=vol[:,i])
  radialprofile[2,i]=np.average(pgas[:,i],weights=vol[:,i])
  radialprofile[3,i]=np.average(Ekin[:,i],weights=vol[:,i])
  radialprofile[4,i]=np.average(vel1[:,i],weights=vol[:,i])
  radialprofile[5,i]=np.average(vel1[:,i],weights=vol[:,i]*rho[:,i])
  radialprofile[6,i]=np.average(vel3[:,i],weights=vol[:,i])
  radialprofile[7,i]=np.average(vel3[:,i],weights=(vol[:,i]*rho[:,i]))
  radialprofile[8,i]=np.average(B1[:,i],weights=vol[:,i])
  radialprofile[9,i]=np.average(B2[:,i],weights=vol[:,i])
  radialprofile[10,i]=np.average(B3[:,i],weights=vol[:,i])
  radialprofile[11,i]=np.average(PB[:,i],weights=vol[:,i])
  radialprofile[12,i]=np.average(Er[:,i],weights=vol[:,i])
  radialprofile[13,i]=np.average(Fr1[:,i],weights=vol[:,i])
  radialprofile[14,i]=np.average(Maxwell[:,i],weights=vol[:,i])
  radialprofile[15,i]=np.average(Reynolds[:,i],weights=vol[:,i])
  radialprofile[16,i]=np.average(kappa[:,i],weights=vol[:,i])
  radialprofile[17,i]=np.average(radacc[:,i],weights=vol[:,i])
  massflux=np.sum(rhov[:,i]*x1area[:,i])/Medd
  massfluxin=np.sum(rhovin[:,i]*x1area[:,i])/Medd
  massfluxout=np.sum(rhovout[:,i]*x1area[:,i])/Medd
  radialprofile[18,i]=massflux
  radialprofile[19,i]=massfluxin
  radialprofile[20,i]=massfluxout

ntlim1=ntheta/4
ntlim2=ntheta*3/4

for i in range(0,nr):
  radialprofile2[1,i]=np.average(rho[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[2,i]=np.average(pgas[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[3,i]=np.average(Ekin[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[4,i]=np.average(vel1[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[5,i]=np.average(vel1[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i]*rho[ntlim1:ntlim2,i])
  radialprofile2[6,i]=np.average(vel3[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[7,i]=np.average(vel3[ntlim1:ntlim2,i],weights=(vol[ntlim1:ntlim2,i]*rho[ntlim1:ntlim2,i]))
  radialprofile2[8,i]=np.average(B1[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[9,i]=np.average(B2[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[10,i]=np.average(B3[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[11,i]=np.average(PB[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[12,i]=np.average(Er[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[13,i]=np.average(Fr1[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[14,i]=np.average(Maxwell[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[15,i]=np.average(Reynolds[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[16,i]=np.average(kappa[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  radialprofile2[17,i]=np.average(radacc[ntlim1:ntlim2,i],weights=vol[ntlim1:ntlim2,i])
  massflux2=np.sum(rhov[ntlim1:ntlim2,i]*x1area[ntlim1:ntlim2,i])/Medd
  massfluxin2=np.sum(rhovin[ntlim1:ntlim2,i]*x1area[ntlim1:ntlim2,i])/Medd
  massfluxout2=np.sum(rhovout[ntlim1:ntlim2,i]*x1area[ntlim1:ntlim2,i])/Medd
  radialprofile2[18,i]=massflux2
  radialprofile2[19,i]=massfluxin2
  radialprofile2[20,i]=massfluxout2



col_names='[1]r    [2]rho    [3]Pg    [4]Ek    [5]Vr    [6]Vr_m    [7]Vphi    [8]Vphi_m    [9]B1    [10]B2    [11]B3    [12]PB    [13]Er    [14]Fr1    [15]Max    [16]Rey    [17]kappa    [18]radacc    [19]Mdot    [20]Mdotin    [21]Mdotout'


first_file='{:05d}'.format(ni)
end_file='{:05d}'.format(no)

np.savetxt('timeave_rprof_'+first_file+'_'+end_file+'.txt',np.transpose(radialprofile),fmt='%4.3e',header=col_names)

np.savetxt('timeave_rprof2_'+first_file+'_'+end_file+'.txt',np.transpose(radialprofile2),fmt='%4.3e',header=col_names)


