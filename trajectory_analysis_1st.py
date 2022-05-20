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
  # --------------------------------------
  dt     = 0.1   # s  
  m      = 10.0  # kg
  v      = 15e3 # m/s
  alt    = 200e3  # m
  gam    = -11   # deg
  theta  = 0     # deg
  Cd = 1.0                 # -
  Cl = 0.0                 # -
  S  = math.pi * 0.4**2    # m^2

  atm  = "atmospheremodel.txt"
  outf = "output.dat"

  # const
  RE = 6.378e6             # m
  g  = 9.8                 # m/s^2

  # --------------------------------------
  # init setting
  # --------------------------------------
  atm_para = np.loadtxt(atm, skiprows=21, unpack=True) # 0:alt, 4:den
  atm_alt = atm_para[0] * 10**3 # m
  atm_den = atm_para[4] * 10**3 # kg/m^3

  gam   = gam   *math.pi/180  # rad
  theta = theta *math.pi/180  # rad
  
  t = 0                # s
  r = alt + RE        # m
  alt = (r - RE)       # m
  rho = interpo(alt, atm_den, atm_alt)     # kg/m^3 
  numQ = 4             # number of q
  Res = np.zeros(numQ) # residual 
  Q   = np.zeros(numQ) # v[m/s], gamma[rad], r[m], theta[rad]
  Q[0] = v
  Q[1] = gam
  Q[2] = r
  Q[3] = theta

  with open(outf,"w") as f:
    f.write("time[s], v[m/s], gam[deg], alt[m], theta[deg], density[kg/m^3]")
    f.write("\n")
  #end
  output(outf, Q, RE, t, rho)

  # --------------------------------------
  # loop
  # --------------------------------------
  while r > RE :
    t  += dt

    Res[0] = - 0.5 * Cd * S * rho * v**2 / m - g * RE**2 / r**2 * np.sin(gam)
    Res[1] =   0.5 * Cl * S * rho / m - (g/v * RE**2 / r**2 - v / r) * np.cos(gam)
    Res[2] = v * np.sin(gam)
    Res[3] = v / r * np.cos(gam)
    
    # 1st-order Euler
    Q[:] = Q[:] + dt * Res[:]
    
    #update
    v     = Q[0]
    gam   = Q[1]
    r     = Q[2]
    theta = Q[3]
    alt = (r - RE)                           # m
    rho = interpo(alt, atm_den, atm_alt)     # kg/m^3 
        
    # output
    output(outf, Q, RE, t, rho)

    # print
    print("{:.2f}km".format(alt/1000))
  #end
#end

def interpo(alt, atm_den, atm_alt):
  rho = 0
  n = len(atm_den)-1

  for i in range(n):
    if atm_alt[i] <= alt and alt <= atm_alt[i+1]:
      a = (atm_den[i+1] - atm_den[i])/(atm_alt[i+1] - atm_alt[i])
      b = atm_den[i] - a * atm_alt[i]
      rho = a * alt + b
      break
    #end
  #end

  if atm_alt[n] < alt:
    rho = atm_den[n]
  elif alt < atm_alt[0] :
    rho = atm_den[0]
  #end

  return rho
#end

def output(fff, Q, RE, t, rho):
  with open(fff,"a") as f:
    f.write("{:.1e}".format(t))
    f.write(" ")
    f.write("{:.6e}".format(Q[0]))
    f.write(" ")
    f.write("{:.6e}".format(Q[1] * 180/math.pi))
    f.write(" ")
    f.write("{:.6e}".format(Q[2] - RE))
    f.write(" ")
    f.write("{:.6e}".format(Q[3] * 180/math.pi))
    f.write(" ")
    f.write("{:.6e}".format(rho))
    f.write(" ")
    f.write("\n")
  #end
#end

# --------------------------------------------------------
# --------------------------------------------------------
if __name__ == "__main__":
    main()