'''
Analyze a lot of files for the calibration of this accelerometer,
provided the data are formatted in a certain specific way.

Uses the optimal parameters found in the last fit as
a guess for the next fit. I find this helps improve the returned optimal params.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('ggplot')

names = 'count, timestamp, temp, x_accel, y_accel, z_accel, x_gyro, y_gyro, z_gyro'

AX_OFFSET = 259.77 # comes from my 7-13 cube offset calculation
GLOC = 9.80258 # local acceleration of gravity from https://www.ngs.noaa.gov/cgi-bin/grav_pdx.prl
dGloc = 0.00002 # if I ever need it

datalist = [] # initialize an empty list of z-sensitivites to print at the end
max_rot_list = [] # similarly, initialize a list of maximum rotational velocities

plot = input('Show plots? Type y for yes: ')

for i in range(10):
    # everything below this indent analyzes one file.
    # as such, path_swing and path_stat need to change based on where the data are
    print('\nTrial', str(i+1)+':\n')
    path_swing = r"C:\Users\mhanr\Desktop\NIST\Data\7-13\W" + str(i+1) + '.txt'
    path_stat = r"C:\Users\mhanr\Desktop\NIST\Data\7-13\T" + str(i+1) + '.txt'

    swing_data = np.genfromtxt(path_swing, names = names, skip_footer = 1)
    stat_data = np.genfromtxt(path_stat, names = names, skip_footer = 1)

    stat_time = stat_data['timestamp']
    stat_time -= stat_time[0] # makes the first entry zero, should help with the fit
    stat_count = stat_data['count']
    stat_ax_o = stat_data['x_accel'] - AX_OFFSET
    stat_gz = stat_data['z_gyro']

    stat_ax_resids = stat_ax_o - (np.average(stat_ax_o))
    stat_gz_resids = stat_gz - (np.average(stat_gz))

    print('Ax summary stats:')
    print('Mean:', np.average(stat_ax_o))
    print('2 Stdev: ', 2*np.std(stat_ax_o), '\n')

    print('Gz summary stats:')
    print('Mean:', np.average(stat_gz))
    print('2 Stdev: ', 2*np.std(stat_gz), '\n')

    if plot == 'y':
        fig, ax = plt.subplots(1, 2)
        ax[0].scatter(stat_time, stat_ax_resids, s = 2, c = 'r')
        ax[0].plot(stat_time, [2*np.std(stat_ax_o)]*len(stat_time), 'k--')
        ax[0].plot(stat_time, [-2*np.std(stat_ax_o)]*len(stat_time), 'k--')
        ax[0].set_title('X Accelerometer')
        ax[0].set_ylabel('Ax - Mean (bits)')
        ax[0].set_xlabel('Time (s)')

        ax[1].scatter(stat_time, stat_gz_resids, s = 2, c = 'b')
        ax[1].plot(stat_time, [2*np.std(stat_gz)]*len(stat_time), 'k--')
        ax[1].plot(stat_time, [-2*np.std(stat_gz)]*len(stat_time), 'k--')
        ax[1].set_title('Z Gyroscope')
        ax[1].set_ylabel('Gz - Mean (bits)')
        ax[1].set_xlabel('Time (s)')
        plt.show()

    gz_offset = np.average(stat_gz)

    swing_time = swing_data['timestamp']
    swing_time -= swing_time[0]
    swing_count = swing_data['count']
    swing_ax_o = swing_data['x_accel'] - AX_OFFSET
    swing_gz_o = swing_data['z_gyro'] - gz_offset

    swing_ax_resids = swing_ax_o - np.average(swing_ax_o)
    swing_gz_resids = swing_gz_o - np.average(swing_gz_o)

    # this is where the fitting begins
    def accel(x, g, k, d, e, f):
        return g * np.exp(-k*x) * np.cos( (d*np.cos(w*x)+e*np.sin(w*x))) + f
    def gyro(t, k, w, a, b, c):
        return np.exp(-k*t) * ( w * ( a * np.cos(w*t) + b * np.sin(w*t) ) ) + c

    # this is where the actual fitting begins
    accel_guesses =  [10696.02,  0.00161577, -0.1102, 0.9386, 8719.41] # you will likely need to change this!
    gyro_guesses = [0.007126082, 5.148, 5753.85, -271.69, 16.24] # and this!

    # sometimes the above guesses really mess up the first one, but they also often help
    # if the first default guess hurts, you can just make an exception if i==0

    popt1, pcov1 = curve_fit(gyro, swing_time, swing_gz_o, gyro_guesses)
    w = popt1[1] # needed for the next fit!
    popt0, pcov0 = curve_fit(accel, swing_time, swing_ax_o, accel_guesses)

    g, ka, d, e, f = popt0
    k, w, a, b, c = popt1

    err0 = dg, dka, dd, de, df = np.sqrt(np.diag(pcov0))
    err1 = dk, dw, da, db, dc = np.sqrt(np.diag(pcov1))

    A = np.average(stat_ax_o)
    dA = np.std(stat_ax_o)

    params = [a,b,d,e,f,g,A]
    uncerts = [da, db, dd, de, df, dg, dA]

    accel_fit = accel(swing_time, *popt0)
    gyro_fit = gyro(swing_time, *popt1)

    if plot == 'y':
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(swing_time, swing_ax_o, 'ro', label = 'x-accel data', markersize = 2)
        ax[0].plot(swing_time, accel_fit, 'k--', label = 'fit')
        ax[0].legend()
        ax[0].set_xlabel('Timestamp (sec)')
        ax[0].set_ylabel("Acceleration (bits)")

        ax[1].plot(swing_time, swing_gz_o, 'go', label = 'ang. accel data', markersize = 2)
        ax[1].plot(swing_time, gyro_fit, 'k--', label = 'fit')
        ax[1].legend()
        ax[1].set_xlabel('Timestamp (sec)')
        ax[1].set_ylabel('Gyroscopic Acceleration (bits)')
        plt.show()

    # here I print summary stats that gaitan does
    Astatic = np.average(stat_ax_o)
    Amin = g*np.cos(np.sqrt(d**2 + e**2))+ f
    theta_0_rad= np.arccos(1-(Astatic-Amin)/(Astatic))
    theta_0 = np.degrees(theta_0_rad)
    Gm = theta_0*w
    Max_Gz = np.sqrt(a**2 + b**2)*w
    Sensitivity = Gm/Max_Gz
    datasheet = 4.4375 # mdps/bit
    difference = ((datasheet / 1000) - Sensitivity) / (datasheet / 1000)

    print('Astatic:', Astatic, 'bits')
    print('Amin:', Amin, 'bits')
    print('Theta 0:', theta_0, 'deg')
    print('Gm = theta0*w:', Gm, 'dps')
    print(Max_Gz, 'bits')
    print(Sensitivity, 'dps/bit')
    print(Sensitivity*1000, 'mdps/bit')
    print('difference: ', "{:.7%}".format(difference))
    print('L = g0/w2 =', GLOC / w**2, 'm')
    plt.show()
    datalist.append(Sensitivity)
    max_rot_list.append(Gm)

# print the data out at the end to copy/paste for analysis
for i in datalist:
    print(i, ',')
print()
for j in max_rot_list:
    print(j, ',')

plot = input('Show more plots? [y for yes]: ')
datasht_sense = 0.004375 #datasheet senstivity
datalist = [i for i in datalist if not np.isnan(i)]

if plot=='y':
    resids = input('Residuals? [y for yes]: ')
    x = np.arange(1, len(datalist)+ 1)

    if resids == 'y':
        upper_tolerance = [np.asarray(np.average(datalist)) *.015]*(len(datalist))
        lower_tolerance = [np.asarray(np.average(datalist)) * (-.015)]*(len(datalist))

        plt.plot(x, upper_tolerance, linestyle = 'dashed', label = 'upper tolerance')
        plt.plot(x, lower_tolerance, linestyle = 'dashed', label = 'lower tolerance')
        plt.plot(x, datalist - np.average(datalist), 'bo', label = 'observed sensitivity')
        plt.plot(x, [0]*(len(x)), '--', )
        plt.title('Residual plot')
        plt.ylabel('Sensitivity (dps/lsb)')
        plt.legend()
        plt.show()

    else:
        expected = [datasht_sense]*len(datalist)
        upper_tolerance = [np.asarray(np.average(datalist)) *1.015]*(len(datalist))
        lower_tolerance = [np.asarray(np.average(datalist)) * (1-.015)]*(len(datalist))  # let the computer do math

        plt.plot(x, expected, linestyle = 'dashed', label = 'datasheet sensitivity')
        plt.plot(x, upper_tolerance, linestyle = 'dashed', label = 'upper tolerance')
        plt.plot(x, lower_tolerance, linestyle = 'dashed', label = 'lower tolerance')
        plt.plot(x, datalist, 'bo', label = 'observed sensitivity')
        plt.plot(x, [np.average(datalist)]*(len(x)), '--', )
        plt.title('Plot')
        plt.ylabel('Sensitivity (dps/lsb)')
        plt.legend()
        plt.show()
        # todo: figure out stdev and plot where nums would be expected vs where they are
        # could be instructive of whether they have systematic error or noise
print(f'mean value is {(np.average(datalist) - datasht_sense)* 100/ datasht_sense :.4g}% off from datasheet sensitivty')

# uncertainty analysis: wip
# import sympy as sp # a python-based symbolic manipulator
# # a list of symbols to represent the independent variables
# symbols = a_s,bs,ds,es,fs,gs,As = sp.symbols('a b d e f g A')
# # the equation for sensitivty that uses our indep variable symbols:
# S_sym = 1/((a_s**2 + bs**2))**(0.5)* sp.acos((gs*sp.cos((ds**2 + es**2)**(0.5)) + fs) / (As))
#
# def grad_func(a,b,d,e,f,g,A):
#     # returns the matrix of numerical partial deriviatves (i.e., the gradient)
#     sym_partial_derivs = [sp.diff(S_sym, var) for var in symbols]
#     grad = []
#     for func in sym_partial_derivs:
#         func = sp.lambdify(symbols, func, 'numpy')
#         grad.append(func(a,b,d,e,f,g,A))
#     return grad
#
# gradient = grad_func(a,b,d,e,f,g,A) # our gradient
# mysum = 0
# for i in range(len(gradient)):
#     mysum += (gradient[i] * uncerts[i])**2
#
# # i wonder if this is actually a better way to calculate it:
# def S(z,y,x,w,v,u,t):
#     # note that this result is in radians per second per bit.
#     return 1/np.sqrt(z**2+y**2)* np.arccos((u*np.cos(np.sqrt(x**2+w**2))+v) / (t))
