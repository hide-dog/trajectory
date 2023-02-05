# --------------------------------------
# import
# --------------------------------------
import numpy as np
import math
import matplotlib.pyplot as pl

# --------------------------------------
# main
# --------------------------------------
fff = "output_1st.dat"
x_py1 = np.loadtxt(fff, skiprows=1, unpack="True")
fff = "output_4th.dat"
x_py4 = np.loadtxt(fff, skiprows=1, unpack="True")
# fff = "result_orbit.dat"
# x_tacode = np.loadtxt(fff, skiprows=2, unpack="True")
# ----------
# pl.plot(x_tacode[1],x_tacode[8],"r-", label="tacode")
pl.plot(x_py1[0],x_py1[1],"g-", label="python_1st")
pl.plot(x_py4[0],x_py4[1],"b-",label="python_4th")
pl.xlabel("Time, s")
pl.ylabel("Velocity, m/s")
pl.grid()
pl.legend(fontsize=14)
pl.savefig("TV.png", format="png", dpi=300)
pl.clf()
pl.close()
# ----------
# pl.plot(x_tacode[1],x_tacode[4],"r-", label="tacode")
pl.plot(x_py1[0],x_py1[3]/1000,"g-", label="python_1st")
pl.plot(x_py4[0],x_py4[3]/1000,"b-", label="python_4th")
pl.xlabel("Time, s")
pl.ylabel("Altitude, km")
pl.grid()
pl.legend(fontsize=14)
pl.savefig("TA.png", format="png", dpi=300)
pl.clf()
pl.close()
# ----------
pl.plot(x_py1[1],x_py1[3]/1000,"g-", label="python_1st")
pl.plot(x_py4[1],x_py4[3]/1000,"b-", label="python_4th")
pl.xlabel("Altitude, km")
pl.ylabel("Velocity, m/s")
pl.grid()
pl.legend(fontsize=14)
pl.savefig("AV.png", format="png", dpi=300)
pl.clf()
pl.close()
# ----------
pl.plot(x_py4[10],x_py4[9],"b-", label="python_4th")
pl.xlabel("Mach number")
pl.ylabel("Reynolds number")
pl.grid()
pl.legend(fontsize=14)
pl.savefig("RM.png", format="png", dpi=300)
pl.clf()
pl.close()
# ----------
pl.plot(x_py4[0],x_py4[12],"r", linestyle="dotted", linewidth = 2.0, label="con")
pl.plot(x_py4[0],x_py4[13],"orange", linestyle="dashed", linewidth = 2.0, label="rad")
pl.plot(x_py4[0],x_py4[14],"k", linestyle="solid", linewidth = 2.0, label="tatal")
pl.xlabel("Time, s")
pl.ylabel("Heat flux, W/m^2")
pl.grid()
pl.legend(fontsize=14)
pl.savefig("Tq.png", format="png", dpi=300)
pl.clf()
pl.close()
# ------------
tsize="14"
fig = pl.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax3 = ax1.twinx()
# ax4 = ax1.twinx()

ax1.plot(x_py4[0],x_py4[3]/1000, color = "k", linestyle="solid", linewidth = 2.0, label="altitude")
ax2.plot(x_py4[0],x_py4[1],      color = "b", linestyle="dotted", linewidth = 2.0, label="velocity")
ax3.plot(x_py4[0],x_py4[11],      color = "g", linestyle="dashed", linewidth = 2.0, label="G")

#-------------------
ax1.grid(True)
ax1.set_xlabel('Tims, s',fontsize=tsize)
ax1.set_ylabel('Altitude, km',fontsize=tsize)
ax2.set_ylabel('Velosity, m/s',fontsize=tsize)
ax3.set_ylabel('Impact acceleration, G',fontsize=tsize)

ax3.spines["right"].set_position(("axes", 1.3))
# ax4.spines["right"].set_position(("axes", 2.4))

#-------------------
# ax1.set_xlim([0.0, 120])
# ax1.set_ylim([0, 200])
# ax2.set_ylim([0, 60])
# ax3.set_ylim([0, 9000])
# ax1.set_yticks(np.arange(36, 44+1, step=1))

""" scale name size  """
ax1.tick_params(labelsize=tsize)
ax2.tick_params(labelsize=tsize)
ax3.tick_params(labelsize=tsize)

fig.legend(fontsize=14, bbox_to_anchor=(0.15, -0.05, 0.5, 1))

""" format adjustment """
pl.tight_layout()


# save
pl.savefig("TAVI.png", format="png", dpi=300)
# reset
pl.clf()
pl.close()

# ----------
RE = 6.378e6   # m          
r = x_py4[3] + RE
theta = (x_py4[4]+90)/180*math.pi
x = np.zeros(len(theta))
y = np.zeros(len(theta))
for i in range(len(theta)):
    x[i] = r[i] * math.cos(theta[i])
    y[i] = r[i] * math.sin(theta[i])

pl.plot(x, y,"b", linestyle="solid", linewidth = 2.0, label="trajectory")

theta = np.linspace(0,360,360) / 180*math.pi
x = np.zeros(len(theta))
y = np.zeros(len(theta))
for i in range(len(theta)):
    x[i] = RE * math.cos(theta[i])
    y[i] = RE * math.sin(theta[i])

pl.plot(x, y,"k", linestyle="solid", linewidth = 2.0, label="Earth")

x *= (RE+100e3)/RE
y *= (RE+100e3)/RE
pl.plot(x, y,"k", linestyle="dashed", linewidth = 0.5, label="Alt 100 km")
x *= (RE+50e3)/(RE+100e3)
y *= (RE+50e3)/(RE+100e3)
pl.plot(x, y,"k", linestyle="dashed", linewidth = 1.0, label="Alt 50 km")

pl.xlabel("x, m")
pl.ylabel("y, m")
pl.xlim([-1.0e6, 0.25e6])
pl.ylim([6.2e6, 6.6e6])
pl.grid()
pl.legend(fontsize=14)
pl.savefig("posi.png", format="png", dpi=300)
pl.clf()
pl.close()