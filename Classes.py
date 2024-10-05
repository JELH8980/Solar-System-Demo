# -----------------------------------------------------------------------------
# Script: Class Definitions for Orbital Simulation
# Author: Ludwig Horvath
# Email: ludhor@kth.se
# Date: 2024-10-05
# Description: This script contains the class definitions for the orbital 
#              simulation, including the Body, System, Orbit, and frame classes. 
#              Each class handles specific functionalities related to celestial 
#              body attributes, orbital mechanics, and graphical representations.
# -----------------------------------------------------------------------------

import numpy as np
from numpy import pi, deg2rad
from numpy import cos as C, sin as S, arctan2, sqrt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import pandas as pd

G = 6.67430*10**(-11)       # Universal gravitational constant

AU = 149597871*10**3 # [m/1 AU]

errtol = 0.01

t0 = 0

sec_per_tunit =  3600*24    # time unit is [day]

dt = sec_per_tunit        # time step [1 day]

earthyear = 365

close = False


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)

class frame():
    def __init__(self, title):
        self.title = title
        self.e_x = None
        self.e_y= None
        self.e_z = None

        self.W, self.i, self.w = 0, 0, 0

        self.T = np.eye(3)

    def R_3(self, x):
        return np.array([[C(x), S(x), 0], [-S(x), C(x), 0], [0, 0, 1]])

    def R_2(self, x):
        return np.array([[C(x), 0, -S(x)], [0, 1, 0], [S(x), 0, C(x)]])

    def R_1(self, x):
        return np.array([[1,0, 0], [0, C(x), S(x)], [0, -S(x), C(x)]])

    def RX_x(self, W, i, w):  
        return self.R_3(w) @ self.R_1(i) @ self.R_3(W)

    def Rx_X(self, W, i, w): 
        return self.RX_x(W,i,w).T

    def orient(self, W, i, w):
        self.W, self.i, self.w = W, i, w
        self.T = self.Rx_X(self.W, self.i, self.w)
        self.e_x = self.T @ self.e_x
        self.e_y = self.T @ self.e_y
        self.e_z = self.T @ self.e_z

class Body:
    def __init__(self, name, mass, radius):
        self.mass = mass       # [kg]
        self.name = name    
        self.radius = radius     # [m]
        self.position = None

class System:
    def __init__(self, name, filepath='Orbital_Elements.csv'):
        self.name = name
        self.filepath = filepath
        self.orbits = []
        self.planets = {}
        self.no_bodies = 0
        self.centerbody = None
        self.colordictionary = {'Sun': '#ffd700',
                                'Mercury': '#b0b0b0',
                                'Venus': '#e1ce7a',
                                'Earth': '#1f77b4',
                                'Mars': '#ff5733',
                                'Jupiter': '#ffad33',
                                'Saturn': '#eedc82',
                                'Uranus': '#7fdbff',
                                'Neptune': '#5073b8',
                                'Pluto': '#be93c5'}

    def create(self):
        OE_data = pd.read_csv(self.filepath, index_col=0)

        planetary_data = pd.read_csv('Planetary_Data.csv', index_col=0)

        planet_names = OE_data.columns

        sun_mass = planetary_data.loc['Mass (10^24 kg)', 'Sun']*10**24

        sun_radius = planetary_data.loc['Radius (km)', 'Sun']*1000

        sun = Body('Sun', sun_mass, sun_radius)

        self.centerbody = sun

        for planet_name in planet_names:
            mass = planetary_data.loc['Mass (10^24 kg)', planet_name]*10**24    # Mass
            radius = planetary_data.loc['Radius (km)', planet_name]*10**3       # radius

            planet_instance = Body(planet_name, mass, radius)

            self.planets[planet_name]= planet_instance

            e =            OE_data.loc['e', planet_name]  # [-] eccentricity of orbit
            a =            OE_data.loc['a', planet_name]  # [au] semi-major axis of orbit
            w =            deg2rad(OE_data.loc['w', planet_name]) # [rad] argument of perihelion of orbit (rel. ecliptive plane) (Small omega w)
            i =            deg2rad(OE_data.loc['i', planet_name]) # [rad] inclination of orbit (rel. ecliptive plane)
            W =            deg2rad(OE_data.loc['W', planet_name]) # [rad] ecliptic longitude of ascending node        (Big omega W)
            M =            deg2rad(OE_data.loc['M', planet_name]) # [rad] mean anomaly of planet

            orbital_elements = {"e": e, "a": a, "i": i, "W": W, "w": w, "M": M}

            new_orbit = Orbit(planet_name + '-Sun', [self.centerbody, planet_instance], 
                              orbital_elements, color=self.colordictionary[planet_name])

            self.orbits.append(new_orbit)
            self.no_bodies += 1


class Orbit():
    def __init__(self, title, bodies, orbital_elements, color='black', N=100, origin=None, comment=''):
        global G
        global AU
        global errtol
        global dt

        self.comment = comment
        self.name = title
        self.N = N
        
        self.origin = np.array(origin) if origin is not None else np.zeros(3)
        
        self.body1 = bodies[0]      # Primary Body Dominating mass
        self.body2 = bodies[1]      # Secondary Body Smaller mass

        self.e = orbital_elements["e"]      # [-]   eccentricity of orbit
        self.a = orbital_elements["a"]      # [au]  semi-major axis of orbit
        self.w = orbital_elements["w"]      # [rad] argument of perihelion of orbit (rel. ecliptive plane) (Small omega w)
        self.i = orbital_elements["i"]      # [rad] inclination of orbit (rel. ecliptive plane)
        self.W = orbital_elements["W"]      # [rad] ecliptic longitude of ascending node 
        self.M = orbital_elements["M"]      # [rad] mean anomaly
        self.mu = G*(self.body1.mass+self.body2.mass)   # [-] Gravitational constant of 2-body system
        
        self.P = self.a*AU*(1-self.e**2)

        self.h = sqrt(self.P*self.mu)

        self.T_period = 2*pi/sqrt(self.mu)*(self.a*AU)**(3/2) if self.e <= 1 else None   # [s] 

        self.E = 0 # Initial guess of eccentric anomaly [rad]

        err = 100   # Initial error set for Newtons Method to solve for Eccentric Anomaly: E

        while err > errtol:
            E_old = self.E

            f = self.E - self.e*S(self.E) - self.M
            fp = 1 - self.e*C(self.E)

            self.E = self.E - f/fp

            err = np.abs(self.E-E_old)

        self.nu = 2*arctan2(sqrt(1+self.e)*S(self.E/2), sqrt(1-self.e)*C(self.E/2))

        self.Mdot = 2*pi/self.T_period

        self.r_pf_vec = None
        self.r_vec = None
        
        # Graphical Properties===================================================START

        self.C = (self.pf(0) + self.pf(pi))/2  
        self.F = self.origin
        self.AP = self.pf(pi)
        self.PE = self.pf(0)

        self.line_AP = np.vstack([self.PE, self.AP]) 
        self.line_Vert = np.vstack([self.pf(-pi/2), self.pf(pi/2)]) 

        

        self.curve = self.pf(np.linspace(0,2*pi,N))

        self.color = color
        
        self.frame = frame(title)

        self.arrow_margin = 1
        self.arrow_prop_dict = dict(mutation_scale=10, arrowstyle='-|>', color=self.color, shrinkA=0, shrinkB=0)

        sel_mat_ex = np.array([[0,0,0],[0,1,0],[0,0,1]])
        
        sel_mat_ey = np.array([[1,0,0],[0,0,0],[0,0,1]])

        self.origin_correction = np.vstack([self.origin, self.origin]).T

        self.origin_correction_x = sel_mat_ex @ np.vstack([self.origin, self.origin]).T

        self.origin_correction_y = sel_mat_ey @ np.vstack([self.origin, self.origin]).T

        self.origin_correction_z = np.vstack([self.origin, self.origin]).T

        self.frame.e_x = self.origin_correction_x + np.array([[self.PE[0]+self.arrow_margin, self.PE[0]+1+self.arrow_margin], [0, 0], [0, 0]])
        self.frame.e_y = self.origin_correction_y + np.array([[0, 0], [self.pf(pi/2)[1]+self.arrow_margin, self.pf(pi/2)[1]+1+self.arrow_margin], [0, 0]])
        self.frame.e_z = self.origin_correction_z + np.array([[0, 0], [0, 0], [0+self.arrow_margin, 1+self.arrow_margin]])
    
        # Graphical Properties======================================================END

        self.frame.orient(self.W,self.i,self.w)
        self.orient_conic()

    def pf(self, nu, shifted=True):
        points = self.P / (1 + self.e * C(nu)) * np.array([C(nu), S(nu), np.zeros_like(nu)])
        return (points.T + self.origin).T if shifted else points

    def update(self):
        T = self.frame.T
        self.nu = 2*arctan2(sqrt(1+self.e)*S(self.E/2), sqrt(1-self.e)*C(self.E/2))
        self.r_pf_vec = self.pf(self.nu, False)
        self.r_vec = T @ self.r_pf_vec + self.origin


    def time_step(self, track=False):
        M_new = self.Mdot*dt + self.M

        err = 100
        E = self.E

        while err > errtol:
            E_old = E

            f = E-self.e*S(E)-M_new
            fp = 1-self.e*C(E)

            E = E - f/fp

            err = np.abs(E-E_old)
        
        self.E = E
        self.M = M_new

        self.update()
        
        if track:
            print('True anomaly: ' + str(self.nu))
            print('Position (ecliptic): ' + str(self.pos))

        return

    def orient_conic(self):

        T = self.frame.T

        self.C = T @ (self.C-self.origin) + self.origin

        self.AP = T @ (self.AP-self.origin) + self.origin

        self.PE = T @ (self.PE-self.origin) + self.origin

        self.curve = T @ (self.curve-self.origin.reshape(3,1)) + self.origin.reshape(3,1)

        self.line_AP = (T @ (self.line_AP.T-self.origin.reshape(3,1))).T + self.origin.reshape(3,1).T

        self.line_Vert = (T @ (self.line_Vert.T-self.origin.reshape(3,1))).T + self.origin.reshape(3,1).T

        self.frame.e_x = self.PE.reshape(3,1) + T @ np.array([[0,1],[0,0],[0,0]])
        self.frame.e_z[:,1] = self.frame.e_z[:,0] + (self.frame.e_x[:,1]-self.frame.e_x[:,0])/np.linalg.norm(self.frame.e_z[:,1]-self.frame.e_z[:,0])

        self.frame.e_y =  self.line_Vert[1,:].T.reshape(3,1) + T @ np.array([[0,0],[0,1],[0,0]])
        self.frame.e_z[:,1] = self.frame.e_z[:,0] + (self.frame.e_y[:,1]-self.frame.e_y[:,0])/np.linalg.norm(self.frame.e_z[:,1]-self.frame.e_z[:,0])

        self.frame.e_z = self.origin_correction_z +  T @ np.array([[0,0],[0,0],[0,1]])
        self.frame.e_z[:,1] = self.frame.e_z[:,0] + (self.frame.e_z[:,1]-self.frame.e_z[:,0])/np.linalg.norm(self.frame.e_z[:,1]-self.frame.e_z[:,0])

        self.arrow_x = Arrow3D(self.frame.e_x[0], self.frame.e_x[1], self.frame.e_x[2], **self.arrow_prop_dict)
        self.arrow_y = Arrow3D(self.frame.e_y[0], self.frame.e_y[1], self.frame.e_y[2], **self.arrow_prop_dict)
        self.arrow_z = Arrow3D(self.frame.e_z[0], self.frame.e_z[1], self.frame.e_z[2], **self.arrow_prop_dict)

