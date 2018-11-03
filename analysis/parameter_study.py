import sys
sys.path.append('.')

from thermofluids import Atmosphere
from components import Rocket, Motor
from trajectory import TwoDTrajectorySolver
from aerodynamics import DragModel
from plotting import plotCurvesWithMax
import matplotlib.pyplot as plt
import numpy as np
import math
from sensitivity_analysis import sensitivity_analysis_apogee

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
print("Baseline Apogee: {:.0f}m".format(solver.apogee))
# Body
# studyParameters = ['body_diameter', 'nose_length', 'body_length', 'boattail_diameter']
# Fins
studyParameters = ['finSpan', 'finTipChord', 'finRootChord', 'finSweepLength']
start = -0.5
stop = 2
step = 0.10
multipliers = np.arange(start, stop+step, step)
print(multipliers)

counter = 1
bApogee = 0

apogeeSensitivityList = []
stabiltiySensitivityList = []
apogeeList = []
stabilityList = []


for parameter in studyParameters:
    
    print("({}/{}) {}".format(counter, len(studyParameters), parameter))
    counter += 1

    # start from baseline rocket every time
    myRocket = Rocket(rocketFilePath)
    myMotor = Motor(motorFilePath)
    
    solver = TwoDTrajectorySolver(
        myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
        dragModel=DragModel.CoefficientsBox
    )

    analysisResults = sensitivity_analysis_apogee(parameter, solver, multipliers)
    
    apogeeSensitivityList.append(analysisResults["apogee_sensitivity"])
    stabiltiySensitivityList.append(analysisResults["stability_sensitivity"])
    apogeeList.append(analysisResults["apogee_list"])
    stabilityList.append(analysisResults["stability_list"])

f, axarr = plt.subplots(2, sharex=True)

for curve in stabilityList:
    axarr[0].plot(multipliers, curve, '-o')
axarr[0].set_title("stability vs rocket parameters")
axarr[0].legend(studyParameters)
axarr[0].set(xlabel="% change in parameter", ylabel="% change in stability")
axarr[0].grid()

for curve in apogeeList:
    axarr[1].plot(multipliers, curve, '-o')
axarr[1].set_title("stability sensitivity vs rocket parameters")
axarr[1].legend(studyParameters)
axarr[1].set(xlabel="% change in parameter", ylabel="% change in apogee")
axarr[1].grid()

plt.show()