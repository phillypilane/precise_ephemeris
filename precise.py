#apg4001
#philly pilane
#ass 3

from math import *
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D

#Constants
GM = 3.986008e14 #product of gravitational constant and mass of earth
OMEGA_DOT_E = 7.292115167e-5 #rotation rate of earth

#Read all lines from broadcast ephermiris file
lines = open('brdc2000.15n','r').readlines()
lines2 = open('igs18540.sp3','r').readlines()

all_lines = [] #empty lists
my_epochs = []

for i,j in enumerate(lines): #enumerated to allow loop to store 8 line segments after sat no.
    stripped = lines[i][:-1]
    all_lines += [stripped]
    if all_lines[i][0:2] == '13':
        my_sat_obs = []
        for k in range(0,8):
            my_sat_obs += [[lines[i+k][:-1].replace("D","e")]]
    else:
        continue
    my_epochs += [my_sat_obs]


radial_pred = []
#extracting values from lines as per gps nav format file structure
for i in range(2,3):
    sat_no = int(my_epochs[i][0][0][0:2])
    TOC = my_epochs[i][0][0][3:22]
    time = my_epochs[i][0][0][12:17].replace(" ",":")
    SV_ClockBias = float(my_epochs[i][0][0][22:41])
    SV_ClockDrift = float(my_epochs[i][0][0][41:60])
    SV_ClockDriftRate = float(my_epochs[i][0][0][60:80])

    IODE = float(my_epochs[i][1][0][3:22]) 
    Crs = float(my_epochs[i][1][0][22:41]) #(meters)
    Delta_n = float(my_epochs[i][1][0][41:60]) #(radians/sec)
    M0 = float(my_epochs[i][1][0][60:80]) #(radians)

    Cuc = float(my_epochs[i][2][0][3:22]) #(radians)
    e_Ecc = float(my_epochs[i][2][0][22:41])
    Cus = float(my_epochs[i][2][0][41:60]) #(radians)
    sqrt_A = float(my_epochs[i][2][0][60:80]) #(sqrt(m))

    Toe = float(my_epochs[i][3][0][3:22]) #sec of week
    Cic = float(my_epochs[i][3][0][22:41]) #(radians)
    OMEGA = float(my_epochs[i][3][0][41:60])#(radians)
    Cis = float(my_epochs[i][3][0][60:80])#(radians)

    i0 = float(my_epochs[i][4][0][3:22]) #(radians)
    Crc = float(my_epochs[i][4][0][22:41]) #(radians)
    small_omega = float(my_epochs[i][4][0][41:60]) #(radians)
    OMEGA_DOT = float(my_epochs[i][4][0][60:80]) #(radians/sec)

    IDOT = float(my_epochs[i][5][0][3:22]) #(radians/sec)
    CODES_L2 = float(my_epochs[i][5][0][22:41])
    GPS_WEEK = float(my_epochs[i][5][0][41:60])
    L2_P = float(my_epochs[i][5][0][60:80])

    SV_accuracy = float(my_epochs[i][6][0][3:22]) #(meters)
    SV_health = float(my_epochs[i][6][0][22:41]) #(bits)
    TGD = float(my_epochs[i][6][0][41:60]) #(seconds)
    IODC = float(my_epochs[i][6][0][60:80])


    time_4plt = []
    xplt=[] #to be plotted
    yplt=[]
    zplt=[]
    for j in range(0,96):
        A = sqrt_A**2 #semi major axis
        n0 = sqrt(GM/(A*A*A)) #computed mean motion
        tk = 900*j - Toe #time from ephemeris reference epoch
        tk_4plt = 900*j
        if tk > 302400:
            tk = tk-604800
        elif tk < -302400:
            tk = tk +604800
        else:
            pass
        n = n0 + Delta_n #corrected mean motion
        Mk = M0 + n*tk #mean anomaly
        Ek = Mk + e_Ecc*sin(Mk)
        iterate = 5
        while iterate > 0:
            Ek = Mk + e_Ecc*sin(Ek)
            iterate+=-1
            
        Qa = sqrt(1-e_Ecc**2)*sin(Ek)
        Qb = 1-e_Ecc*cos(Ek)
        Qc = cos(Ek)-e_Ecc
        Qd = 1-e_Ecc*cos(Ek)
        vk = atan((Qa/Qb)/(Qc/Qd)) #true anomaly

        phi_k = vk + small_omega #argument of latitude

        #second harmonic perturbations
        duk = Cus*sin(2*phi_k)+Cuc*cos(2*phi_k) #argument of latitude corr
        drk = Crc*cos(2*phi_k)+Crs*sin(2*phi_k) #radius corr
        dik = Cic*cos(2*phi_k)+Cis*sin(2*phi_k) #corr to inclination
        
        uk = phi_k +duk #corrected argument of latitude
        rk = A*(1 -e_Ecc*cos(Ek))+drk #corrected radius
        ik = i0 + dik + (IDOT)*tk #corred inclination

        xk_ = rk*cos(uk) #x position in orbital plane
        yk_ = rk*sin(uk) #y position in orbital plane

        omega_k = OMEGA + ((OMEGA_DOT-OMEGA_DOT_E)*tk)-(OMEGA_DOT_E*Toe)


        xk = xk_*cos(omega_k)-yk_*cos(ik)*sin(omega_k)
        yk = xk_*sin(omega_k)+yk_*cos(ik)*cos(omega_k)
        zk = yk_*sin(ik)

        xplt.append(xk)
        yplt.append(yk)
        zplt.append(zk)
        

        radial_calc = sqrt(((0-xk)**2)+((0-yk)**2)+((0-zk)**2))
        radial_pred += [radial_calc]
        time_4plt += [tk_4plt]

###############################  SP3  ##################################

    

all_lines2 = []
radial_prec = []

xplot = []
yplot = []
zplot = []

for i in range(len(lines2)):
    stripped2 = lines2[i][:-1]
    all_lines2 += [stripped2]
    if all_lines2[i][0:4] == 'PG13':
        xxx = float(all_lines2[i][5:18])*1000
        yyy = float(all_lines2[i][19:32])*1000
        zzz = float(all_lines2[i][33:46])*1000

        xplot.append(xxx)
        yplot.append(yyy)
        zplot.append(zzz)
    else:
        continue

    radial_calc2 = sqrt(((0-xxx)**2)+((0-yyy)**2)+((0-zzz)**2))
    radial_prec += [radial_calc2]
        
diff_list = []
for i in range (len(radial_prec)):
    diff = radial_pred[i] - radial_prec[i]
    diff_list += [diff]

plt.plot(time_4plt, diff_list, 'b') #graph


small_diff = []
time4_small = []

for i in range(len(diff_list)):
    test = diff_list[i]
    test_time = time_4plt[i]
    if -2.0 <= test <= 2.0:
        small_diff.append(test)
        time4_small.append(test_time)
    else:
        continue

####################################### RMS ################################
numerator = []
for i in range(len(small_diff)):
    big_sig = (small_diff[i] - np.mean(small_diff))**2
    numerator.append(big_sig)

RMS = sqrt(sum(numerator)/len(small_diff))

print ('The RMS for the interval with differences less than 2m is: ', round(RMS,4),'m')

plt.plot(time4_small, small_diff, 'r', linewidth = 2.5)
plt.xlabel('Time (seconds)')
plt.ylabel('Radial Vector Differences (meters)')
plt.title('Difference between Broadcast Ephermris and Precise Ephermeris over time')#graph
plt.show()
 
fig=pylab.figure()
ax=Axes3D(fig)
ax.scatter(xplot,yplot,zplot,c=u'r')
ax.scatter(xplt,yplt,zplt,c=u'g')
pyplot.show()

