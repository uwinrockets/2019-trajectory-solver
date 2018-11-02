import math

class Atmosphere:
  
    sea_pressure = 101.325 * math.pow(10, 3) # Pa
    sea_temperature = 288.15 # K
    lapse_rate = 0.0065 # K/m
    grav_accel = 9.80665 # m/s/s
    gas_constant_universal = 8.31447 # J/mol/K
    molar_mass = 0.0289644 # kg/mol
    standard_gravity = 9.80665 #m/s/s
    earth_mean_radius =  6371 * math.pow(10, 3) #m
    gamma = 1.4
    gas_constant_air = 287 #J/kg/K
    # for sutherland
    c1 = 1.458 * math.pow(10, -6)
    s = 110.4
  
    def __init__(self, altitude):
        self.altitude = altitude
        self.density = None
        self.dyn_viscosity = None
        self.kin_viscosity = None
        self.gravity = None
        self.temperature = self.sea_temperature
        self.pressure = self.sea_pressure
        self.speedOfSound = None

        self.updateAltitude(altitude)
    
    def updateAltitude(self, altitude):
        self.altitude = altitude
        # these must be updated first
        self.temperature = self.__calculateTemperature()
        self.pressure = self.__calculatePressure()
        self.density = self.__calculateDensity()
        
        self.dyn_viscosity = self.__calculateDynViscosity()
        self.kin_viscosity = self.dyn_viscosity/self.density
        self.gravity = self.__calculateGravity()
        self.speedOfSound = self.__calculateSpeedOfSound()

    # https://en.wikipedia.org/wiki/Density_of_air#Altitude 
    def __calculateTemperature(self):
        return self.sea_temperature - (self.lapse_rate * self.altitude)

    # https://en.wikipedia.org/wiki/Density_of_air#Altitude  
    def __calculatePressure(self):
        return self.sea_pressure * math.pow(\
            (1 - (self.lapse_rate * self.altitude)/self.sea_temperature),\
            (self.grav_accel * self.molar_mass)/(self.gas_constant_universal * self.lapse_rate)\
        )

    # https://en.wikipedia.org/wiki/Density_of_air#Altitude
    def __calculateDensity(self):
        return (self.pressure * self.molar_mass)/(self.gas_constant_universal * self.temperature)

    # https://www.cfd-online.com/Wiki/Sutherland%27s_law
    def __calculateDynViscosity(self):
        return (self.c1 * math.pow(self.temperature, 3.0/2.0)) / (self.temperature + self.s)

    def __calculateGravity(self):
        return self.standard_gravity * math.pow(self.earth_mean_radius/(self.earth_mean_radius + self.altitude), 2)
    
    def __calculateSpeedOfSound(self):
        return math.sqrt(self.gamma * self.gas_constant_air * self.temperature)