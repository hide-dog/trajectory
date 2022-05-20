
# ---------------------------------------------------------------------------------------------------
# memo
# Rd is not accurate
#
# ---------------------------------------------------------------------------------------------------

# --------------------------------------
# import
# --------------------------------------
import numpy as np
import math

# --------------------------------------
# main
# --------------------------------------
def main():
  # --------------------------------------
  # setting
  # ---------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------
  dt    = 0.1   # s  
  m     = 10.0  # kg
  v     = 15e3  # m/s
  alt   = 200e3 # m
  gam   = -11   # deg
  theta = 0     # deg

  Cd = 1.0                 # -
  Cl = 0.0                 # -
  S  = math.pi * 0.4**2    # m^2
  l  = 0.8                 # m
  
  atm  = "atmospheremodel.txt"
  outf = "output.dat"
  timeInt = 1              # 1:1st-euler, 2:4RK
  # ---------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------
  # const
  RE = 6.378e6             # m
  g  = 9.8                 # m/s^2
  R  = 8.314               # J/K/mol
  M = np.array([0.028, 0.014, 0.032, 0.016])     # kg/mol, N2, N, O2, O
  Cv = np.array([2.5, 1.5, 2.5, 1.5])     # kg/mol, N2, N, O2, O
  Cp = np.array([3.5, 2.5, 3.5, 2.5])     # kg/mol, N2, N, O2, O
  # ---------------------------------------------------------------------------------------------------

  # --------------------------------------
  # init setting
  # ---------------------------------------------------------------------------------------------------
  atm_para = np.loadtxt(atm, skiprows=21, unpack=True) # 0:alt, 1:O, 2:N2, 3:O2, 4:den, 5:tem, 6:N
  atm_alt = atm_para[0] * 10**3 # m
  atm_den = atm_para[4] * 10**3 # kg/m^3
  atm_tem = atm_para[5]         # K
  atm_rat = np.array([atm_para[2], atm_para[6], atm_para[3], atm_para[1]])   # -, N2, N, O2, O

  gam   = gam   *math.pi/180  # rad
  theta = theta *math.pi/180  # rad
  
  chem = np.zeros(4)     # volume ratio of Chemical species
  
  t = 0                # s
  r = alt + RE         # m
  alt = (r - RE)       # m
  rho, tem, chem = interpolate(alt, atm_den, atm_tem, atm_alt, atm_rat, chem)     # kg/m^3, K

  numQ = 4                 # number of q
  Res = np.zeros(numQ)     # residual 
  k4  = np.zeros((4,numQ)) # k4
  Q   = np.zeros(numQ)     # v[m/s], gamma[rad], r[m], theta[rad]
  Q[0] = v
  Q[1] = gam
  Q[2] = r
  Q[3] = theta
  
  phyval = np.zeros(11) # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma
  input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv)

  with open(outf,"w") as f:
    # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma
    f.write("time[s], v[m/s], gam[deg], alt[m], theta[deg], density[kg/m^3], temperature[K], ")
    f.write("pressure[Pa], dynamic pressure[Pa], Re, Ma")
    f.write("\n")
  #end
  output(outf, phyval)
  # ---------------------------------------------------------------------------------------------------

  # --------------------------------------
  # loop
  # ---------------------------------------------------------------------------------------------------
  while r > RE :
    # time step
    t  += dt

    if timeInt == 1:
      # 1st-order Euler
      Q = first_euler(numQ, Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)
    elif timeInt == 2:
      # 4th-order Runge-Kutta
      Q = forth_runge_kutta(numQ, Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g, k4)
    else:
      print(" --------------------- ")
      print(" change timeInt !")
      print(" meow !")
      print(" --------------------- ")
      exit()
    #end
    
    # update
    r     = Q[2]
    alt = (r - RE)                                         # m
    rho, tem, chem = interpolate(alt, atm_den, atm_tem, atm_alt, atm_rat, chem)     # kg/m^3, K

    # output
    input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv)
    output(outf, phyval)

    # print
    print("{:.2f} km".format(alt/1000))
  #end
  # ---------------------------------------------------------------------------------------------------
#end

# --------------------------------------
# 1st-order euler
# --------------------------------------
def first_euler(numQ, Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g):
  Res = cal_res(Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)
  
  Q[:] = Q[:] + dt * Res[:]
  
  return Q
#end

# --------------------------------------
# 4th-order Runge-Kutta
# --------------------------------------
def forth_runge_kutta(numQ, Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g, k4):
  
  iRes = np.copy(Res)
  
  k4[0][:] = cal_res(Q, iRes, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)
  iRes[:] = Res[:] + 0.5*dt * k4[0][:]
  
  k4[1][:] = cal_res(Q, iRes, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)
  iRes[:] = Res[:] + 0.5*dt * k4[1][:]

  k4[2][:] = cal_res(Q, iRes, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)
  iRes[:] = Res[:] +     dt * k4[2][:]

  k4[3][:] = cal_res(Q, iRes, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g)

  Q[:] = Q[:] + dt * 1/6 *(k4[0][:] + 2*k4[1][:] + 2*k4[2][:] + k4[3][:])
  
  return Q
#end

# --------------------------------------
# calculate residual
# --------------------------------------
def cal_res(Q, Res, dt, RE, atm_den, atm_tem, atm_alt, atm_rat, chem, Cd, Cl, S, m, g):
  v     = Q[0]
  gam   = Q[1]
  r     = Q[2]
  alt = (r - RE)                           # m
  rho, tem, chem = interpolate(alt, atm_den, atm_tem, atm_alt, atm_rat, chem)     # kg/m^3, K

  Res[0] = - 0.5 * Cd * S * rho * v**2 / m - g * RE**2 / r**2 * np.sin(gam)
  Res[1] =   0.5 * Cl * S * rho / m - (g/v * RE**2 / r**2 - v / r) * np.cos(gam)
  Res[2] = v * np.sin(gam)
  Res[3] = v / r * np.cos(gam)
  return Res
#end

# --------------------------------------
# interpolate
# --------------------------------------
def interpolate(alt, atm_den, atm_tem, atm_alt, atm_rat, chem):
  rho = 0.0
  tem = 0.0
  n = len(atm_den)-1

  for i in range(n):
    if atm_alt[i] <= alt and alt <= atm_alt[i+1]:
      rho = interpolate_xxx(alt, atm_alt, atm_den, i)
      tem = interpolate_xxx(alt, atm_alt, atm_tem, i)
      chem[0] = interpolate_xxx(alt, atm_alt, atm_rat[0], i)
      chem[1] = interpolate_xxx(alt, atm_alt, atm_rat[1], i)
      chem[2] = interpolate_xxx(alt, atm_alt, atm_rat[2], i)
      chem[3] = interpolate_xxx(alt, atm_alt, atm_rat[3], i)
      chem[:] /= np.sum(chem) 
      break
    #end
  #end

  if atm_alt[n] < alt:
    rho = atm_den[n]
    tem = atm_tem[n]
    chem[:] = atm_rat.T[n]
  elif alt < atm_alt[0] :
    rho = atm_den[0]
    tem = atm_tem[0]
    chem[:] = atm_rat.T[0]
  #end

  return rho, tem, chem
#end

def interpolate_xxx(alt, atm_alt, atm_xxx, i):
  a = (atm_xxx[i+1] - atm_xxx[i])/(atm_alt[i+1] - atm_alt[i])
  b = atm_xxx[i] - a * atm_alt[i]
  temp = a * alt + b
  return temp
#end
# --------------------------------------
# calculate physical value
# --------------------------------------
def input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv):
  # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma
  v     = Q[0]
  gam   = Q[1] * 180/math.pi
  theta = Q[3] * 180/math.pi

  Rd = R / (np.dot(chem[:], M[:]))
  k  = np.dot(chem[:], Cp[:] / Cv[:])

  a = (k * Rd * tem)**0.5
  mu = 1.82e-5 * (tem/293.15)**1.5 * (293.15 + 117)/(tem + 117)

  phyval[0]  = t
  phyval[1]  = v
  phyval[2]  = gam
  phyval[3]  = alt
  phyval[4]  = theta
  phyval[5]  = rho
  phyval[6]  = tem
  phyval[7]  = rho * Rd * tem
  phyval[8]  = 0.5 * rho * v**2
  phyval[9]  = rho * v * l / mu
  phyval[10] = v / a

  return phyval
#end

# --------------------------------------
# output
# --------------------------------------
def output(fff, val):
  with open(fff,"a") as f:
    for i in range(len(val)):
      f.write("{:.6e}".format(val[i]))
      f.write(" ")
    #end
    f.write("\n")
  #end
#end

# --------------------------------------------------------
# --------------------------------------------------------
if __name__ == "__main__":
    main()