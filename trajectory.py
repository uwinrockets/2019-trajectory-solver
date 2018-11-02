from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import integrate
from components import Rocket, Motor
from aerodynamics import getDragCoefficient, DragModel

'''
constant mass, c_d, gravity
'''
def f(t, z, params):
    x, vx, y, vy = z

    rocket, motor, atmosphere, dragModel, launchAngle, launchrailLength, constantDrag = params

    mass = rocket.geometry["motorless_mass"] + motor.getMass(t)

    # updates the atmospheric properties
    atmosphere.updateAltitude(y)
    density = atmosphere.density
    g = atmosphere.gravity
    a = atmosphere.speedOfSound
    
    theta = None
    if (y <= launchrailLength):
        theta = launchAngle
    else:
        theta = math.atan(vy/vx)

    V = math.sqrt(vx**2+vy**2)
    mach = V/a

    c_d = getDragCoefficient(rocket, motor.isPowered(t), atmosphere, mach, dragModel, constantDrag)

    drag = (1.0/2.0)*density*(V**2)*rocket.geometry["referenceArea"]*c_d['all']
    thrust = motor.getThrust(t)
    
    derivatives = [
        vx,
        -1.0*(1.0/mass)*drag*np.cos(theta) + (1.0/mass)*thrust*np.cos(theta),
        vy,
        -1.0*(1.0/mass)*drag*np.sin(theta) + (1.0/mass)*thrust*np.sin(theta) - g,
    ]
  
    return derivatives

# event to stop the solver when the ground is hit
def hit_ground(t, y): return y[2]

class TwoDTrajectorySolver:
    def __init__(self, rocket, motor, atmosphere, initialConditions, launchAngle, launchrailLength,\
        hitGround=True, runtime=200, dragModel=DragModel.TacticleMissileDesign, constantDrag = 0.3):

        self.rocket = rocket
        self.motor = motor
        self.atmosphere = atmosphere
        self.initialConditions = initialConditions
        self.launchAngle = launchAngle
        self.launchrailLength = launchrailLength
        self.hitGround = hitGround
        self.runtime = runtime
        self.dragModel = dragModel
        self.constantDrag = constantDrag

        self.parameters = [self.rocket, self.motor, self.atmosphere, self.dragModel,\
                            self.launchAngle, self.launchrailLength, self.constantDrag]

        hit_ground.terminal = self.hitGround
        hit_ground.direction = -1
    
    def solve(self):
        self.soln = integrate.solve_ivp(fun=lambda t, y: f(t, y, self.parameters),\
            t_span=[0., self.runtime], y0=self.initialConditions, events=hit_ground, max_step=0.1)
        
        self.v = self.soln["y"][0]
        self.vx = self.soln["y"][1]
        self.y = self.soln["y"][2]
        self.vy = self.soln["y"][3]
        self.t = self.soln["t"]
        self.apogee = np.max(self.y)
        self.apogee_time = self.t[np.argmax(self.y)]

    def post_process(self):
        self.machVsTime = np.zeros(len(self.t))
        self.dragForceVsTime = np.zeros(len(self.t))
        self.stagTempVsTime = np.zeros(len(self.t))
        self.thrustVsTime = np.zeros(len(self.t))
                
        # dict of lists containing the different drag coefficients vs time
        # can be called like dragCoefficientsVsTime['all'] or dragCoefficientsVsTime['fin_friction']
        self.dragCoefficientsVsTime = {}
        firstLoop = True

        # i is the times at which the solution has been computed
        for i in range(len(self.t)):
            self.atmosphere.updateAltitude(self.y[i])

            V = math.sqrt(self.vx[i]**2 + self.vy[i]**2)
            self.machVsTime[i] = V / self.atmosphere.speedOfSound
            # print("temp {:.2f} mach {:.2f} altitude {:.2f}".format(self.atmosphere.temperature, self.machVsTime[i], self.y[i]))
            self.stagTempVsTime[i] = self.atmosphere.temperature * (1 + ((1.4-1)/2)*self.machVsTime[i]**2)
            self.thrustVsTime[i] = self.motor.getThrust(i)

            dragCoefficients = getDragCoefficient(
                self.rocket, self.motor.isPowered(self.t[i]),\
                self.atmosphere, self.machVsTime[i], self.dragModel,\
                constantDrag=self.constantDrag
            )

            if firstLoop:
                firstLoop = False
                for key in dragCoefficients:
                    self.dragCoefficientsVsTime[key] = np.zeros(len(self.t))
            
            for key in dragCoefficients:
                self.dragCoefficientsVsTime[key][i] = dragCoefficients[key]
            
            self.dragForceVsTime[i] = (1./2.) * self.atmosphere.density * V**2 * dragCoefficients['all'] * self.rocket.geometry["referenceArea"]

        self.accelY = np.diff(self.vy)/np.diff(self.t)
        self.accelX = np.diff(self.vx)/np.diff(self.t)
            
        self.maxMach = np.max(self.machVsTime)
        self.maxTemp = np.max(self.stagTempVsTime)
        self.maxDrag = np.max(self.dragForceVsTime)
        self.maxThrust = np.max(self.thrustVsTime)
    
    def printPostProcess(self):
        print("Apogee: {:.2f}m at {:.2f}s".format(self.apogee, self.apogee_time))
        print("Max Mach: {:.2f}".format(self.maxMach))
        print("Max Temp: {:.2f}K at Mach {:.2f} and {:.0f}m".format(
            self.maxTemp, self.machVsTime[np.argmax(self.stagTempVsTime)], self.y[np.argmax(self.stagTempVsTime)]))
        print("Max Drag: {:.2f}N".format(self.maxDrag))
        print("Max AccelY: {:.2f}m/s/s".format(np.max(self.accelY)))
        print("Max AccelX: {:.2f}m/s/s".format(np.max(self.accelX)))
        print("Max Thrust: {:.2f}N".format(self.maxThrust))