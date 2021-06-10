
'''Analyze the readings from the accelerometer/gyroscope to fit
gyroscopic acceleration in the z-direction and linear acceleration in the x-direction.'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path = input('file name: ')

# turn all the data into numpy arrays
names = 'count, timestamp, x_accel, y_accel, z_accel, x_gyro, y_gyro, z_gyro'
all_data = np.genfromtxt(path, names = names, skip_footer = 1)

# this calibration factor came from the 6/8 calibration of device taped to wall,
# where I took one reading with the x-axis facing up and one
# with the z-axis facing down
M_S_2_PER_BIT = 0.0005994292667111093
x_sensitivity = 16351.4205 # g per lsb in for the x-axis

# so for linear acceleration, we can convert the bits to m/s^2
timestamp = all_data['timestamp']
x_accel = all_data['x_accel']/x_sensitivity
y_accel = all_data['y_accel']/x_sensitivity
z_accel = all_data['z_accel']/x_sensitivity
z_gyro_bits = all_data['z_gyro']

# define functions designed for scipy's curve fit that can guess paramters
def x_accel_to_fit(t, g, ka, d, e, f):
    return g *np.exp(-ka*t)* np.cos( (d*np.cos(w*t)+e*np.sin(w*t)) ) + f
def theta_accel(t, k, w, a, b, c):
    return np.exp(-k*t) * ( w * ( a * np.cos(w*t) + b * np.sin(w*t) ) ) + c

# this is where the actual fitting begins
popt0_guesses =  [1, .1, 1, 1, np.average(x_accel)]
popt1_guesses = [.01, 4, 2500, 2500, np.average(z_gyro_bits)]

popt1, pcov1 = curve_fit(theta_accel, timestamp, z_gyro_bits, popt1_guesses)

w = popt1[1] # needed for the next fit!
popt0, pcov0 = curve_fit(x_accel_to_fit, timestamp, x_accel)

# plot the fits
accel_fit = x_accel_to_fit(timestamp, *popt0)
rad_accel_fit = theta_accel(timestamp, *popt1)
fig, ax = plt.subplots(2,1, figsize = (6,12))
ax[0].plot(timestamp, x_accel, c = 'dodgerblue', label = 'x-accel data')
ax[0].plot(timestamp, accel_fit, c= 'orangered', label = 'fit')
ax[0].legend()
ax[0].set_xlabel('Timestamp (sec)')
ax[0].set_ylabel('Acceleration (m$/s^{2}$)')

ax[1].plot(timestamp, z_gyro_bits, c = 'dodgerblue', label = 'ang. accel data')
ax[1].plot(timestamp, rad_accel_fit, c= 'orangered', label = 'fit')
ax[1].legend()
ax[1].set_xlabel('Timestamp (sec)')
ax[1].set_ylabel('Gyroscopic Acceleration (bits)')

'''
--------------------------------------------------------------------------------
The rest of this script prints the fitting parameters and some calculations.
Comment it all out except for plt.show() if you don't care about numbers.
--------------------------------------------------------------------------------
'''
# print and guess the fit parameters:
print('x acceleration:')
titles = ['x accel amplitude', 'damping coef', 'a', 'b', 'height']
units = ["g's", '/s', 'bits', 'bits', "g's"]
for i in range(len(titles)):
    print(titles[i], ': ', popt0[i] , '+/-', pcov0[i][i]**0.5, units[i])
print('\nrad accel:')
titles2 = ['damping coef', 'ang freq', 'd', 'e', 'height']
units2 = ['/s', 'rad/s', 'bits', 'bits', 'bits']
for i in range(len(titles2)):
    print(titles2[i], ': ', popt1[i] , '+/-', pcov1[i][i]**0.5, units2[i])

# this cell is the back-of-envelope calculations Dr. Gaitan showed in our teams meeting 6/7
L = g_m_s2/popt1[1]**2
max_gz_bits = popt1[1]*np.sqrt(popt1[2]**2+popt1[3]**2)
Axmin = min(x_accel)*x_sensitivity
theta0 = np.cos(Axmin/x_sensitivity)
gz_max_rad_s = theta0*popt1[1]
max_gz_bits = popt1[1]*np.sqrt(popt1[2]**2+popt1[3]**2)
z_sensitivity = np.degrees(gz_max_rad_s)/ (max_gz_bits)
datasheet_sensitivity = .004375

print('L: ', L, 'm')
print('Max Gz:', max_gz_bits, 'bits')
print('Ax min:', Axmin)
print('theta0:', theta0, 'rad,', np.degrees(theta0), 'deg')
print('Gz max:', gz_max_rad_s, 'rad/s, ', np.degrees(gz_max_rad_s), 'deg/s')
print('z-sensitivity: ', z_sensitivity, 'dps/bits')
print('datasheet: ', datasheet_sensitivity)
print('ratio: ', (datasheet_sensitivity - z_sensitivity)/ (z_sensitivity))

plt.show() # don't comment me!
