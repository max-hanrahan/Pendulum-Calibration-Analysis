#!/usr/bin/env python
# coding: utf-8

'''
A python version of auto_analysis.ipynb. You should be able to type in the name of the
file to be analyzed and your guess of several parameters, and it should output its guess for the natural angular frequency.
So far it only attempts to fit gyroscopic acceleration in the z-direction and linear acceleration in the x-direction
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# TODO: add damping term for x-acceleration
# todo: convert bit reading into quantities with proper units
# todo: add more fits
path = input('file name (type "def" to use default): ')
if path == 'def':
    path = 'AccelGyro Test.txt'

# turn all the data into numpy arrays
# there is likely a much faster way to do this but I don't care (yet)
all_data = np.genfromtxt(path, delimiter = '\t', skip_footer = 1)

# for some reason count = all_data[:,0] reads it as a list of nan's
# so I'm just using arange to avoid that
count = np.arange(0,1000, 1)
timestamp = all_data[:,1]
temp = all_data[:,2]
x_accel_bits = all_data[:,3]
y_accel_bits = all_data[:,4]
z_accel_bits = all_data[:,5]
x_gyro_bits = all_data[:,6]
y_gyro_bits = all_data[:,7]
z_gyro_bits = all_data[:,8]

# constants:
g = 9.80665 # from https://physics.nist.gov/cgi-bin/cuu/Value?gn

# these data were taken with the device stationary atop the workbench
# we therefore have a guess calibration factor
M_S_2_PER_BIT = g/np.average(z_accel_bits)

# so for linear acceleration, we can convert the bits to m/s^2
x_accel = x_accel_bits*M_S_2_PER_BIT
y_accel = y_accel_bits*M_S_2_PER_BIT
z_accel = z_accel_bits*M_S_2_PER_BIT

'''
The hardest part of the fitting is the guessing.
You might have to run it many times, manually adjusting these parameters every time.
'''

### parameter guesses
THETA_0 = np.radians(float(input('initial angle guess (deg): ')))
k_guess = float(input('k-factor guess (for linear damping): '))
frequency = float(input('guess angular frequency: '))

# define functions designed for scipy's curve fit that can guess omega

def x_accel_to_fit(t, frq, theta_init):
    # optimizes fit for linear acceleration in the x-direction by changing theta and omega
    return (g*np.cos(theta_init * np.cos(frq * t)))*M_S_2_PER_BIT + np.average(x_accel) # add mean for better fit

def theta_accel(t, omega, k_factor):
    # optimizes the fit gyroscopic acceleration by changing omega and k
    return -THETA_0 * omega**2 * np.cos(omega*t) * np.exp(-k_factor * t) + np.average(z_gyro_bits) # add mean for better fit


# this is where the actual fitting begins
parameters0 = [frequency, THETA_0]
parameters1 = [frequency, k_guess]

# popt is the list of optimal parameters, as determined by curve_fit
# pcov is the covariance matrix, whose diagonal entries are the variance of each parameter
popt0, pcov0 = curve_fit(x_accel_to_fit, timestamp, x_accel, parameters0)
popt1, pcov1 = curve_fit(theta_accel, timestamp, z_gyro_bits, parameters1)

accel_fit = x_accel_to_fit(timestamp, *parameters0)
rad_accel_fit = theta_accel(timestamp, *parameters1)
fig, ax = plt.subplots(2,1, figsize = (6,12))
ax[0].plot(timestamp, x_accel, c = 'dodgerblue', label = 'x-accel data')
ax[0].plot(timestamp, accel_fit, c= 'orangered', label = 'fit')
ax[0].legend()
ax[0].set_xlabel('Timestamp (sec)')
ax[0].set_ylabel('Acceleration (bits)')

ax[1].plot(timestamp, z_gyro_bits, c = 'dodgerblue', label = 'ang. accel data')
ax[1].plot(timestamp, rad_accel_fit, c= 'orangered', label = 'fit')
ax[1].legend()
ax[1].set_xlabel('Timestamp (sec)')
ax[1].set_ylabel('Gyroscopic Acceleration (bits)')


'''Keep re-running with changed parameters until the fit no longer looks like trash.
Our ultimate goal is to find $\omega$, which we find like so:'''

# omega found from fitting the x-acceleration, with its std dev as uncertainty:
print(popt0[0], '+/-', pcov0[0][0]**(1/2))

#omega found from fitting the angular acceleration, with its std dev as uncrtainty:
print(popt1[0], '+/-', pcov1[0][0]**(1/2))
plt.show()
