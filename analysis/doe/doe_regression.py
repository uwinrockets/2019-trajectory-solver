from __future__ import division

import sys
sys.path.append('.')

from components import Rocket

import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn import linear_model
import itertools
from sklearn.metrics import mean_squared_error, r2_score
from mpl_toolkits.mplot3d import Axes3D

from scipy import optimize

def getRegressionXs(doe_coeff):
    counter = 0
    x_list = []
    xs = np.zeros(doe_coeff.shape)

    for row in doe_coeff:
        for i in range(len(parameterKeys)):
            xs[counter][i] = baselineParameters[parameterKeys[i]] * row[i]
            
        regression_row = []
        
        # x_i
        for i in range(nParameters):
            regression_row.append(xs[counter][i])
        
        # x_i**2
        for i in range(nParameters):
            regression_row.append(xs[counter][i]**2)
        
        # x_i*x_j
        for i in range(0, nParameters):
            for j in range(1, nParameters-i):
                regression_row.append(xs[counter][i]*xs[counter][j])
        
        x_list.append(regression_row)
        
        counter += 1
    return x_list

rocketFilePath = "2019_rocket_doe.csv"
myRocket = Rocket(rocketFilePath)

parameterKeys = ['finRootChord', 'finTipChord', 'finSpan', 'finSweepLength']
baselineParameters = {}
for key in parameterKeys:
    baselineParameters[key] = myRocket.geometry[key]

doe_coeff = np.loadtxt('analysis/doe/lhc_n4_1000samples.csv')
y_apogee = np.loadtxt('analysis/doe/apogee_list_lhs1000.csv')
y_stability = np.loadtxt('analysis/doe/stability_list_lhs1000.csv')

nSamples = len(y_apogee)
nParameters = len(parameterKeys)
xs_regression = getRegressionXs(doe_coeff)
apogeeReg = linear_model.LinearRegression()
stabilityReg = linear_model.LinearRegression()

apogeeReg.fit(xs_regression, y_apogee)
stabilityReg.fit(xs_regression, y_stability)

# check the accuracy of the linear model #############
doe_coeff_test = np.loadtxt('analysis/doe/lhc_n4_100samples.csv')
y_apogee_test = np.loadtxt('analysis/doe/apogee_list_lhs100_test.csv')
y_stability_test = np.loadtxt('analysis/doe/stability_list_lhs100_test.csv')
xs_test = getRegressionXs(doe_coeff_test)
y_apogee_pred = apogeeReg.predict(xs_test)
y_stability_pred = stabilityReg.predict(xs_test)

print('Apogee Variance score: %.2f' % r2_score(y_apogee_test, y_apogee_pred))
print('Stability Variance score: %.2f' % r2_score(y_stability_test, y_stability_pred))
# ####################################################

# http://reliawiki.org/index.php/Response_Surface_Methods_for_Optimization
class DesirabilityFunction:
    def __init__(self, target, weight, lowerLimit=None, upperLimit=None):
        if (lowerLimit!=None) and (upperLimit==None):
            self.mode = 'maximize'
        elif (lowerLimit==None) and (upperLimit!=None):
            self.mode = 'minimize'
        elif (lowerLimit!=None) and (upperLimit!=None):
            self.mode = 'range'
            
        self.T = target
        self.w = weight
        self.L = lowerLimit
        self.U = upperLimit
        
    def eval(self, y):
        if self.mode == 'maximize':
            if y<self.L:
                return 0
            elif y>=self.L and y<=self.T:
                return math.pow((y-self.L)/(self.T-self.L), self.w)
            elif y>self.T:
                return 1
        elif self.mode == 'minimize':
            if y<self.T:
                return 1
            elif y>=self.T and y<=self.U:
                return math.pow((self.U-y)/(self.U-self.T), self.w)
            elif y>self.U:
                return 0
        elif self.mode == 'range':
            if y<self.L:
                return 0
            elif y>=self.L and y<=self.T:
                return math.pow((y-self.L)/(self.T-self.L), self.w)
            elif y>=self.T and y<=self.U:
                return math.pow((self.U-y)/(self.U-self.T), self.w)
            elif y>self.U:
                return 0
  
# assuming they are all of equal importance
def calculateDesirability (response, d_list, importance):
    product = 1
    for (d, r, i) in itertools.izip(d_list, response, importance):
        product *= math.pow(d.eval(r), i)
    return math.pow(product, 1.0/len(d_list))

d_apogee = DesirabilityFunction(9100, 1, lowerLimit=5000)
d_stability = DesirabilityFunction(2.5, 1, lowerLimit=0, upperLimit=5)


# Optimization #######################################
def optimizeWrapper(x, d_list):
    reg_coeff = getRegressionXs(np.array([x]))
    apogee_response = apogeeReg.predict(reg_coeff)
    stability_response = stabilityReg.predict(reg_coeff)
    
    return 1 - calculateDesirability([apogee_response, stability_response], d_list, [1, 1])
    
x0 = (
    0,
    0,
    0,
    0
)

class Attemps:
    def __init__(self):
        self.attemps = []
    def add(self, item):
        self.attemps.append(item)
   
class MyBounds(object):
    def __init__(self, xmax=[1,1,1,1], xmin=[0,0,0,0] ):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        test1 = x[0] > x[1] # root > tip
        test2 = x[0] > x[2] # root > span
        test3 = 1.2217304764 > math.atan(x[3]/x[2]) # sweep angle under 70deg
        return tmax and tmin and test1 and test2 and test3
        
def print_fun(x, f, accepted):
    if f < 1. and accepted:
        # myAttemps.add(x)
        print("f(x): {:.2f} x: {}".format(f, x))
   
func = lambda x: optimizeWrapper(x, [d_apogee, d_stability])
minimizer_kwargs = {"method": "BFGS"}
# minimizer_kwargs = {"method": "Nelder-Mead"}
mybounds = MyBounds()
# myAttemps = Attemps()
ret = optimize.basinhopping(func, x0, minimizer_kwargs=minimizer_kwargs,
                   niter=5000, stepsize=0.05, accept_test=mybounds, callback=print_fun, disp=True)

print(ret)
print("global minimum: x = {}, f(x0) = {:.2f}" .format(ret.x, ret.fun))
reg_coeff = getRegressionXs(np.array([ret.x]))
apogee_response = apogeeReg.predict(reg_coeff)
stability_response = stabilityReg.predict(reg_coeff)
print("predicted apogee: {}".format(apogee_response))
print("predicted stability: {}".format(stability_response))
