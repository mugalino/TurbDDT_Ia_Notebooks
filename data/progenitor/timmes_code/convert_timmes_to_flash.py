import argparse
import numpy as np
from scipy import interpolate

parser = argparse.ArgumentParser(description="Convert Timmes code's progenitors into FLASH format.")
parser.add_argument('name', type=str,
                    help='path to the input file')
args = parser.parse_args()

fp_old = args.name
fp_new = args.name+"_flash"

# Format Reminder:
# Timmes: "i     mass_grav       radius      pressure    density     temp        enuc        ener        egrav"
# FLASH: "# radius         dens           temp        c12          ne22". The second line just contains the amount of rows as integer.

data_old = np.loadtxt(fp_old,skiprows=1)

# Stepping in Timmes is varying a lot based on the integrator. Convert to uniform stepping as in FLASH fiducual progenitor via linear interpolation.
r = data_old[:,2]
r_interp = np.arange(4e5,r.max(),8e5)
nsteps = r_interp.shape[0]
ylist = [] # each entry corresponds to a column (to be) interpolated
ylist.append(r_interp) # r_interp is by definition already interpolated
for i in range(3,data_old.shape[1]):
    y = data_old[:,i]
    f = interpolate.interp1d(r, y)
    ylist.append(f(r_interp))
ylist.append(0.5*np.ones(nsteps)) # we just assume 50/50 C/O
ylist.append(np.zeros(nsteps)) # we assume no metallicity


# Honor FLASH fluff
## 1.) Impose floors
rhofluff = 1e-3
tfluff = 1e3#3e7
ylist[2] = np.where(ylist[2]<rhofluff,rhofluff,ylist[2])
ylist[3] = np.where(ylist[3]<tfluff,tfluff,ylist[3])
## 2.) Impose a temperature decay at rho=1e6 to fluff temperature. This visually reproduces the behaviour of the fiducial model.
idx = np.abs(ylist[2]-1e6).argmin() 
grad = (tfluff-ylist[3][idx])/(r_interp[-1]-r_interp[idx])
ylist[3][idx:] = ylist[3][idx]+grad*(r_interp[idx:]-r_interp[idx])




with open(fp_new, 'w') as f:
    f.write('# radius     dens           temp        c12          ne22\n')
    f.write('%i'%nsteps+'\n')
    for i in range(nsteps):
        for j in [0,2,3,-2,-1]: # indices to map to FLASH format [we have some offset 3 from Timmes format in interpolation step, see above]
            f.write('%.5E'%ylist[j][i]+'  ')
        if i<nsteps-1:
            f.write('\n')
