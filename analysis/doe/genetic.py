from __future__ import division

import sys
sys.path.append('.')

from components import Rocket

import array
import random
import copy
import numpy as np
import math
from sklearn import linear_model
import itertools

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import matplotlib.pyplot as plt
import networkx

# Creating RSM ###########################

# second order fit polynomial
# input:
# (x1, x2, ..., xn)
# output:
# (x1, x2, ..., xn, x1**2, x2**2, ..., xn**2, x1*x2, x1*xn, ..., x2*xn)
def getRegressionXs(xs):
    counter = 0
    x_list = []

    nParameters = len(xs)
    
    regression_row = []
    
    # x_i
    for i in range(nParameters):
        regression_row.append(xs[i])
    
    # x_i**2
    for i in range(nParameters):
        regression_row.append(xs[i]**2)
    
    # x_i*x_j
    for i in range(0, nParameters):
        for j in range(1, nParameters-i):
            regression_row.append(xs[i]*xs[j])
    
    
    return regression_row

rocketFilePath = "2019_rocket_doe.csv"
myRocket = Rocket(rocketFilePath)

parameterKeys = ['finRootChord', 'finTipChord', 'finSpan', 'finSweepLength']
baselineParameters = {}
for key in parameterKeys:
    baselineParameters[key] = myRocket.geometry[key]

doe_coeff = np.loadtxt('analysis/doe/lhc_n4_1000samples.csv')
y_apogee = np.loadtxt('analysis/doe/apogee_list_lhs1000.csv')

nSamples = len(y_apogee)
nParameters = len(parameterKeys)
xs_regression = []
for row in doe_coeff:
    xs_regression.append(getRegressionXs(row))
apogeeReg = linear_model.LinearRegression()

apogeeReg.fit(xs_regression, y_apogee)

##########################################

# Desirability Functions #################
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


d_apogee = DesirabilityFunction(9100, 1, lowerLimit=0, upperLimit=9500)
d_stability = DesirabilityFunction(2.5, 1, lowerLimit=-5, upperLimit=5)
         
##########################################

creator.create("FitnessMulti", base.Fitness, weights=(1.0, 1.2))
creator.create("Individual", np.ndarray, fitness=creator.FitnessMulti)

def func(ind):
    x_regression = getRegressionXs(ind)
    x_regression = np.array(x_regression).reshape(1, -1)
    currentApogee = apogeeReg.predict(x_regression)
    
    # sketchy but leave this for now
    # its because the doe variable were scaled down by 20 inches
    x_stab = []
    for x in ind:
        x_stab.append(x*20/39.37)

    # rocket = Rocket(rocketFilePath)
    myRocket.geometry['finRootChord'] = x_stab[0]
    myRocket.geometry['finTipChord'] = x_stab[1]
    myRocket.geometry['finSpan'] = x_stab[2]
    myRocket.geometry['finSweepLength'] = x_stab[3]
    myRocket.updateGeometry()
    currentStability = myRocket.staticStability

    # print("ind: {} apogee: {} stability: {}".format(ind, currentApogee, currentStability))
    
    return (d_apogee.eval(currentApogee), d_stability.eval(currentStability))

def feasible(ind):
    """Feasability function for the individual. Returns True if feasible False
    otherwise."""
    # implement geometry tests here
    
    test1 = ind[0] > ind[1] # root > tip
    test2 = ind[0] > ind[2] # root > span
    test3 = 1.2217304764 > math.atan(ind[3]/ind[2]) # sweep angle under 70deg
    
    if test1 and test2 and test3:
        return True
    return False

# Structure initializers

# number of parameters
IND_SIZE=4
# number of desirability factors
NDIM = 2
BOUND_LOW, BOUND_UP = 0.0, 1.0

toolbox = base.Toolbox()

toolbox.register("attr_float", random.uniform, BOUND_LOW, BOUND_UP)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_float, n=IND_SIZE)
toolbox.register("evaluate", func)
toolbox.decorate("evaluate", tools.DeltaPenality(feasible, 0))
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0)
toolbox.register("mutate", tools.mutPolynomialBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0, indpb=1.0/IND_SIZE)
toolbox.register("select", tools.selNSGA2)

stats = tools.Statistics()
stats.register("pop", copy.deepcopy)

# history = tools.History()
# toolbox.decorate("mate", history.decorator)
# toolbox.decorate("mutate", history.decorator)

toolbox.pop_size = 100
toolbox.max_gen = 1000
toolbox.mut_prob = 0.2

def run_ea(toolbox, stats=None, verbose=False):
    pop = toolbox.population(n=toolbox.pop_size)
    hof = tools.ParetoFront(similar=np.array_equal)
    # history.update(pop)

    if stats != None:
        stats.register("avg", np.mean, axis=0)
        # stats.register("std", np.std, axis=0)
        # stats.register("min", np.min, axis=0)
         #stats.register("max", np.max, axis=0)
        
    return algorithms.eaMuPlusLambda(pop, toolbox, mu=toolbox.pop_size, 
                                     lambda_=toolbox.pop_size, 
                                     cxpb=1-toolbox.mut_prob,
                                     mutpb=toolbox.mut_prob, 
                                     stats=stats, 
                                     ngen=toolbox.max_gen, 
                                     verbose=verbose,
                                     halloffame=hof)

def main():
    # random.seed(10)
    # ind1 = toolbox.individual()
    # ind1.fitness.values = func(ind1)
    # print(ind1)

    random.seed(21)

    results, logbook = run_ea(toolbox, stats=stats, verbose=True)

    par = []
    d1_list = []
    d2_list = []

    for gen in logbook:
        ind = gen['avg']
        
        d1, d2 = func(ind)
        d1_list.append(d1)
        d2_list.append(d2)

        if feasible(ind):
            par.append(ind)

    count = np.linspace(1, len(logbook), len(logbook))

    plt.plot(count, d1_list, label="apogee")
    plt.plot(count, d2_list, label="stability")
    plt.show()

    # Is not meaningful when the number of generations is so high
    # graph = networkx.DiGraph(history.genealogy_tree)
    # graph = graph.reverse()     # Make the grah top-down
    # colors = [toolbox.evaluate(history.genealogy_history[i])[0] for i in graph]
    # networkx.draw(graph, node_color=colors)
    # plt.show()
    
    print(logbook[-1]['avg'])

if __name__ == "__main__":
    main()