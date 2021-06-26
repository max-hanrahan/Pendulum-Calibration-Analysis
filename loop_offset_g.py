'''
Converts bits to m/s^2 using two text files: one with device stationary facing
down and one with the device stationary facing up.

Loops through 10 similarly named-files!

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('ggplot')

path_to_csv = input("What do you want to name the csv? (just don't type .csv)") + '.csv'

plot = input("Show plots? Type y/n:")

# write the headers exactly once
with open(path_to_csv, 'w') as file_object:
    towrite = 'Sensitivity' +','+ 'Offset' +','+ 'Error' +','+ 'Temperature' +','+ '\n'
    file_object.write(towrite)

for i in range(16):
    # sorry, I know machine dependencies are bad. todo: fix this someday
    path_up = r"C:\Users\mhanr\Desktop\NIST\Data\6-25\Cu"+ str(i+1) + '.txt'
    path_down = r"C:\Users\mhanr\Desktop\NIST\Data\6-25\Cd"+ str(i+1) + '.txt'

    # read in the data
    names = 'count, timestamp, temp, x_accel, y_accel, z_accel, x_gyro, y_gyro, z_gyro'
    all_data_up = np.genfromtxt(path_up, names = names, skip_footer = 1)
    all_data_down = np.genfromtxt(path_down, names = names, skip_footer = 1)

    count = np.arange(0, 1000, 1)

    timestamp_up = all_data_up['timestamp']
    timestamp_down = all_data_down['timestamp']

    x_up_accel_bits = all_data_up['x_accel']
    x_down_accel_bits = all_data_down['x_accel']
    y_up_accel_bits = all_data_up['y_accel']

    y_down_accel_bits = all_data_down['y_accel']
    z_up_accel_bits = all_data_up['z_accel']
    z_down_accel_bits = all_data_down['z_accel']
    # we'll ignore gyroscopic acceleration for now!

    g = 9.80258 # from https://www.ngs.noaa.gov/cgi-bin/grav_pdx.prl using my local coordinates and altitude
    delta_g = 0.00002 # may use this later

    x_sensitivity = (np.average(x_up_accel_bits) - np.average(x_down_accel_bits))/2
    x_offset = (np.average(x_down_accel_bits) + np.average(x_up_accel_bits))/2
    M_S_2_PER_BIT = g/x_sensitivity
    print(M_S_2_PER_BIT)

    # for linear acceleration, we can convert the bits to m/s^2
    x_accel_up = x_up_accel_bits * M_S_2_PER_BIT
    y_accel_up = y_up_accel_bits * M_S_2_PER_BIT
    z_accel_up = z_up_accel_bits * M_S_2_PER_BIT

    x_accel_down = x_down_accel_bits * M_S_2_PER_BIT
    y_accel_down = y_down_accel_bits * M_S_2_PER_BIT
    z_accel_down = z_down_accel_bits * M_S_2_PER_BIT

    print('x senisitivity: ', x_sensitivity, 'bits/g')
    print('x offset: ', x_offset, 'bits/g')
    x_sigma = np.sqrt((np.std(x_up_accel_bits) / 2)**2 + (np.std(x_down_accel_bits) / 2)**2) # two stdevs added in quadrature
    print('x uncert (k = 1):', x_sigma, ',' ,'{:.2g}'.format(x_sigma / x_sensitivity*100), '%', 'fractional')
    print()
    temp_up = np.average(all_data_up['temp'])
    temp_down = np.average(all_data_down['temp'])
    temp_avg = (temp_up + temp_down) / 2.

    # write the sensitivity, offset
    with open(path_to_csv, "a") as file_object:
        towrite = str(x_sensitivity)+','+ str(x_offset) +','+str(x_sigma) +','+ str(temp_avg)+ ','+ '\n'
        file_object.write(towrite)

    # now make the plots
    if plot == 'y':
        # put mean and 2 sigma on each plot in one for loop:
        accel_list_up = [x_accel_up, y_accel_up, z_accel_up]
        accel_list_down = [x_accel_down, y_accel_down, z_accel_down]

        accel_lists = [accel_list_up, accel_list_down]

        fig, ax = plt.subplots(2,3, figsize = (20, 12))
        ax[0,0].scatter(timestamp_up, x_accel_up, c = 'r', s= 2)
        ax[0,1].scatter(timestamp_up, y_accel_up, c = 'g', s= 2)
        ax[0,2].scatter(timestamp_up, z_accel_up, c = 'b', s= 2)

        ax[1,2].scatter(timestamp_down, z_accel_down, c = 'b', s= 2)
        ax[1,0].scatter(timestamp_down, x_accel_down, c = 'r', s= 2)
        ax[1,1].scatter(timestamp_down, y_accel_down, c = 'g', s= 2)

        for x in ax.flat:
            x.set(xlabel='Timestamp', ylabel='Acceleration ($m/s^2$)')

        for x in range(2):
            for i in range(3):
                textstr = '\n'.join((
                r'mean=%.2f' % (np.average(accel_lists[x][i])),
                r'$2\sigma=%.2f$' % (np.std(2* accel_lists[x][i]), )))
                ax[x,i].text(1, np.max(accel_lists[x][i]), textstr, fontsize=14,
                verticalalignment='top')
        plt.show()

# now that that's finished, we analyze the csv with all the data!
names = 'Sensitivity,Offset,Sigma,Temperature'
data = np.genfromtxt(path_to_csv, skip_header = 1, names = names, delimiter = ',')

import os
cwd = os.getcwd()
print('The data are now in', str(path_to_csv), 'in', str(cwd))
plot = input('Show more plots? (y/n):')
if plot == 'y':
    datalist = data['Sensitivity']
    x = np.arange(1, len(datalist)+ 1)
    upper_tolerance = [np.asarray(np.average(datalist)) *1.015]*(len(datalist))
    lower_tolerance = [np.asarray(np.average(datalist)) * (1-.015)]*(len(datalist))  # let the computer do math

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    
    ax1.plot(x, upper_tolerance, linestyle = 'dashed', label = 'upper tolerance')
    ax1.plot(x, lower_tolerance, linestyle = 'dashed', label = 'lower tolerance')
    ax1.plot(x, datalist, 'bo', label = 'observed sensitivity')
    ax1.plot(x, [np.average(datalist)]*(len(datalist)), '--', )
    ax1.set(ylabel = 'Sensitivity (dps/lsb)', title = 'Sensitivity vs Time')
    ax1.legend()

    ax2.errorbar(x, data['Sensitivity'], yerr = np.std(data['Sensitivity'])/np.sqrt(len(data)),
                 fmt = 'bo', capsize = 4,)
    ax2.set(ylabel = 'Sensitivity (bits/g)', title = 'Sensitivity')

    ax3.errorbar(x, data['Offset'], yerr = np.std(data['Offset'])/np.sqrt(len(data)),
                 fmt = 'ro', capsize = 4, )
    ax3.set(ylabel = 'Offset (bits)', title = 'Offset')

    # plot offset consistency against temp
    ax4.plot(data['Temperature'], data['Offset'], 'o')
    ax4.set(xlabel = 'Temperature (Celsius)', ylabel = 'Offset (bits)', title = 'Offset vs Temperature')

    plt.show()
