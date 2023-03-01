import numpy as np
from scipy.integrate import odeint

# constants 
G = 1
c = 1

class blackHole(object):
    """
    Use ODEs to solve simple orbits for BBH systems.
    Solutions assumed in the frame of the more massive m1.
    
    """
    def __init__(self, mass, distance):
        self.mass = mass
        self.distance = distance
           
    def r_init(self, distance):
        return np.array([0,self.distance,0])

    
# two body eom (circular)
def two_body_eom(state, time, G, m1, m2):
    """
    A function to create a circular orbit for BH2, where BH1 is assumed stationary
    state : initial conditions, or the previous state of the system at time t-dt (setup to use scipy.odeint)
    time : array of times
    G : grav constant, assume  G=1
    m1 : mass of BH whose frame we are in - should be larger than m2
    m2 : mass of orbiting BH
    """
    # assume in frame of obj 1
    # find loc and vel of m2 relative to m1
    r_vec = state[3:6]-state[0:3]
    r_mag = np.linalg.norm(r_vec)
    vel = state[9:12] - state[6:9]

    # calc dvdt and return
    a = -(m1) * G * r_vec / (r_mag)**3

    return np.concatenate((state[6:9], vel, [0,0,0], a))

# two body eom with pn1
def two_body_eom_pn(state, time, G, m1, m2):
    """
    A function to create an orbit with PN0 and PN1 for BH2, where BH1 is assumed stationary
    state : initial conditions, or the previous state of the system at time t-dt (setup to use scipy.odeint)
    time : array of times
    G : grav constant, assume  G=1
    m1 : mass of BH whose frame we are in - should be larger than m2
    m2 : mass of orbiting BH
    """
    
    # constants
    mtot = m1 + m2
    eta = (m1*m2)/mtot**2
    
    # assume in frame of obj 1
    # location of BH2 rel to BH1
    r_vec = state[3:6]-state[0:3]
    r_mag = np.linalg.norm(r_vec)
    n_hat = r_vec/r_mag
    
    # speed of BH2 rel to BH1 (BH1 assumed stationary)
    v_vec = state[9:12] - state[6:9]
    v_mag = np.linalg.norm(v_vec)
    v_hat = v_vec/v_mag
    
    # Find PN1 then calc dvdt with PN0 and PN1
    PN1 = (4+2*eta)*(G*mtot/r_mag)*n_hat-(1+3*eta)*v_mag**2*n_hat+(3/2)*eta*(np.dot(n_hat,v_hat))**2*n_hat+(4-2*eta)*np.dot(n_hat,v_hat)*v_vec
    dvdt = (G*mtot/r_mag**2)*(-n_hat+(1/c**2)*PN1)

    return np.concatenate((state[6:9], v_vec, [0,0,0], dvdt))

def orbit_init(bh1,bh2):
    """
    Get masses, radii, and velocities for each BH object.
    Return IC in format for ode solver
    """
    m1 = bh1.mass
    m2 = bh2.mass
    mtot = m1 + m2
    r1 = bh1.r_init(bh1.distance)
    v1 = np.array([0,0,0])
    r2 = bh2.r_init(bh2.distance)
    v2 = np.array([np.sqrt(mtot/bh2.distance), 0, 0])

    ic = np.concatenate((r1,r2,v1,v2))
    return m1, m2, ic


def orbit_ode(eom, init_cond, t, G, m1, m2):
    """
    Solve 2 body eom input by user using scipy's odeint
    eom : a 2 body eom that returns 3D drdt and dvdt
    init_cond : array with initial r and drdt of BH
    t : array of times to solve ode
    G : grav constant, assume  G=1
    m1 : mass of BH whose frame we are in - should be larger than m2
    m2 : mass of orbiting BH
    """
    ysol = odeint(eom,init_cond,t, args=(G,m1, m2))

    x2 = ysol[:,3]
    y2 = ysol[:,4]
    z2 = ysol[:,5]  

    return x2,y2,z2
