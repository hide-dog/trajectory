
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
  # hayabusa condition
  dt     = 0.1   # s
  m      = 16.3  # kg
  v      = 11.8e3 # m/s
  alt    = 200e3  # m
  gam    = -12.7   # deg
  theta  = 0     # deg
  Cd = 1.147                 # -
  Cl = 0.0                 # -

  l  = 0.4                 # length, m
  S  = math.pi * (l/2)**2    # m^2
  
  Rn = 0.2                 # nose radius, m
  tem_w = 300              # temperature at wall, K
  # ----
  
  atm  = "atmospheremodel.txt" # atmosphere model from ncep
  outf = "output_4th.dat"          # output
  timeInt = 1                  # time int, 1:1st-euler, 2:4RK
  # ---------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------
  # const
  RE = 6.378e6             # m
  g  = 9.81                 # m/s^2
  R  = 8.314               # J/K/mol
  M = np.array([0.028, 0.014, 0.032, 0.016])     # kg/mol, N2, N, O2, O
  Cv = np.array([2.5, 1.5, 2.5, 1.5])     # kg/mol, N2, N, O2, O
  Cp = np.array([3.5, 2.5, 3.5, 2.5])     # kg/mol, N2, N, O2, O
  rhos = 1.230             # density at sea level, kg/m^3
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
  
  phyval = np.zeros(15) # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma, G[G], qc[W/m^2], qr[W/m^2], qt[W/m^2]
  input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv, tem_w, Rn, rhos, Cd, m, g)

  with open(outf,"w") as f:
    # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma
    f.write("time[s], v[m/s], gam[deg], alt[m], theta[deg], density[kg/m^3], temperature[K], ")
    f.write("pressure[Pa], dynamic pressure[Pa], Re, Ma, G[G], qc[MW/m^2], qr[MW/m^2], qt[MW/m^2]")
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
    input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv, tem_w, Rn, rhos, Cd, m, g)
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
def input_phyval(phyval, t, Q, alt, rho, tem, l, chem, M, R, Cp, Cv, tem_w, Rn, rhos, Cd, m, g):
  # time[s], v[m/s], gamma[rad], alt[m], theta[rad], rho[kg/m^3], tem[K], pre[Pa], dynamic pre[Pa], Re, Ma, G[G], qc[W/m^2], qr[W/m^2], qt[W/m^2]
  v     = Q[0]
  gam   = Q[1] * 180/math.pi
  theta = Q[3] * 180/math.pi

  Rd = R / (np.dot(chem[:], M[:]))
  Cp_hat = np.dot(chem[:], Cp[:])
  Cv_hat = np.dot(chem[:], Cv[:])
  k  = Cp_hat / Cv_hat
  
  a = (k * Rd * tem)**0.5
  mu = 1.82e-5 * (tem/293.15)**1.5 * (293.15 + 117)/(tem + 117)

  ig = 0.5 * rho * v**2 * (l**2/4*math.pi) * Cd/ m/ g
  
  #---------------------------
  # N.H. KEMP, F.R. RIDDELL, 
  # "Heat Transfer to Satellite Vehicles Re-entering the Atmosphere",
  # Journal of Jet Propulsion. 27 (1957) 132–137. doi:10.2514/8.12603.
  hs  = 0.5*v**2 + Cp_hat*tem
  hw  = Cp_hat*tem_w
  hw0 = Cp_hat*tem_w
  qc = 1.1035*10**4 / Rn**0.5 * (rho/rhos)**0.5 * (v/7925)**3.15 * (hs-hw)/(hs-hw0) # 100 MW/m^2
  # print(1.1035*10**4 / Rn**0.5)
  # print((rho/rhos)**0.5)
  # print((v/7.925)**3.15)
  # print((hs-hw)/(hs-hw0))
  # print(qc/100)
  #---------------------------
  
  #---------------------------
  # Tauber, Michael E., and Kenneth Sutton.
  # "Stagnation-point radiative heating relations for Earth and Mars entries." 
  # Journal of Spacecraft and Rockets 28.1 (1991): 40-42.
  a_tauber = 1.072*10**6 * v**(-1.88) * rho**(-0.325)
  fEV = fev_tauber(v)
  qr = 4.736*10**4 * Rn**a_tauber * rho**1.22 *fEV # 100 MW/m^2
  #---------------------------

  qt = qc + qr

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
  phyval[11] = ig
  phyval[12] = qc/100
  phyval[13] = qr/100
  phyval[14] = qt/100

  return phyval
#end

# --------------------------------------
# output
# --------------------------------------
def output(fff, val):
  if val[3] > 0 :
    with open(fff,"a") as f:
      for i in range(len(val)):
        f.write("{:.6e}".format(val[i]))
        f.write(" ")
      #end
      f.write("\n")
    #end
  #end
#end

# --------------------------------------
# tauber Radiative heating velocity functions for Earth
# --------------------------------------
def fev_tauber(v):
  fev = 0
  tau_v = [9000,  9250,  9500,  9750, 10000, 10250,
          10500, 10750, 11000, 11500, 12000, 12500,
          13000, 13500, 14000, 14500, 15000, 15500, 16000]
  tau_fEv = [1.5, 4.3, 9.7, 19.5, 35, 55, 
              81, 115, 151, 238, 359, 495,
              660, 850, 1065, 1313, 1550, 1780, 2040]
  n = len(tau_v)
  
  for i in range(n):
    if tau_v[i] <= v and v <= tau_v[i+1]:    
      fev = interpolate_xxx(v, tau_v, tau_fEv, i)
      break
    #end
  #end
  if v < tau_v[0]:
    fev = tau_fEv[0]
  elif v > tau_v[n-1]:
    fev = tau_fEv[n-1]
  #end

  return fev
#end

# --------------------------------------------------------
# --------------------------------------------------------
if __name__ == "__main__":
    main()