import sys
sys.path.append('.')

from thermofluids import Atmosphere
from components import Rocket, Motor
from trajectory import TwoDTrajectorySolver
from aerodynamics import DragModel

from scipy import optimize
from scipy.optimize import Bounds
import numpy as np
import math
import pprint

# additional parameters the solver needs
atmosphere = Atmosphere(0)
launchAngle = 84*math.pi/180.0 #rad
launchrailLength = 8.0 #m

# Initial conditions
# y_0 is 0.01 instead of 0 so that the hit_ground event works. look into doing it another way
# x, vx, y, vy
y0 = [0., 0., 0.1, 0.]

rocketFilePath = "2019_rocket.csv"
myRocket = Rocket(rocketFilePath)

motorFilePath = "Cesaroni_N5800.eng"
myMotor = Motor(motorFilePath)

solver = TwoDTrajectorySolver(
    myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
    dragModel=DragModel.CoefficientsBox
)

solver.solve()
print("Baseline Apogee {:.2f}".format(solver.apogee))
print("Baseline Static Stability {:.2f}".format(myRocket.staticStability))

# Class to record the different combinations attempted
class optimizationResults():
    def __init__(self):
        self.counter = 0
        self.geometry = []
        self.apogeeList = []
    def update(self, x, apogee):
        self.geometry.append(x)
        self.apogeeList.append(apogee)
        self.counter += 1
        
def solverOptimizationWrapper(geometryValues, geometryKeys, results, targetApogee=9144):
    # print(geometryValues)    
    for key, value in zip(geometryKeys, geometryValues):
        solver.rocket.geometry[key] = value
        if value < 0:
            return 99999999999
            
    solver.rocket.geometry["body_length"] = 130/39.37 - solver.rocket.geometry["nose_length"] - solver.rocket.geometry["boattail_length"]
    solver.rocket.updateGeometry()
    
    solver.solve()
    target = targetApogee - solver.apogee
    
    results.update(geometryValues, solver.apogee)
    
    print("Iteration {}. Stability {:2f} Apogee {:.2f}m Minimize {:.2f}".format(results.counter, solver.rocket.staticStability, solver.apogee, target))
    return target

optimizationParameters = [
    "nose_length",        # x0
]

bnds = (
    (6/39.37, 100/39.37),  # nose length
)

initialConditions = (
    30/39.37,             # nose length
)

# Inequality constraints of the form f(x) > 0
cons = (
    # enforce the total length of the rocket to be under 130 inches
    {"type": "ineq", "fun": lambda x: solver.rocket.geometry["totalLength"]*39.37 - 130, "label": "length requirement > 130in"},
    # {"type": "ineq", "fun": lambda x: solver.rocket.staticStability - 2, "label": "stability"},    
)

resultsReader = optimizationResults()

fun = lambda x: solverOptimizationWrapper(x, optimizationParameters, resultsReader)
results = optimize.minimize(
    fun, initialConditions, method="SLSQP", bounds=bnds, constraints=(), options={'ftol': 1e-03, 'maxiter': 100}
)

# Post Processing ##################################

print(results)

# convert back to inches/degrees and print results
print('')
print("Results")

for key, value in zip(optimizationParameters, results["x"]):
    if key == "finLESweep":
        print("{}: {:.2f}deg".format(key, value * 180/math.pi))
    else:
        print("{}: {:.2f}cm".format(key, value * 100))
        print("{}: {:.2f}in".format(key, value * 39.37))

print('')        
print("Verifying constraints...")
print("Center of Pressure: {:.2f}".format(solver.rocket.centerPressure))
print("Center of Gravity: {:.2f}".format(solver.rocket.centerGravity))
print("Static Stability: {:.2f}".format(solver.rocket.staticStability))