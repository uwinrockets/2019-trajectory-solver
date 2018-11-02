from __future__ import division
from components import Motor
import math
from enum import Enum

class DragModel(Enum):
    TacticleMissileDesign = 1
    CoefficientsBox = 2
    Constant = 3
    Interpolate = 4

def compressibilityCorrection(coefficient, mach):
    if mach >= 1.1:
        return coefficient / math.sqrt(mach**2-1)
    elif mach <= 0.8:
        return coefficient / math.sqrt(1 - mach**2)
    elif 0.8 < mach < 1.1:
        return coefficient / math.sqrt(1 - 0.8**2)

class CoefficientsTacticleMissile:
    def __init__(self):
        pass

    @staticmethod
    def c_do_body(mach, speedOfSound, density, isPowered,\
        body_length, nose_length, boattail_length,\
        body_diameter, nose_diameter, motor_diameter):    
        
        # meters to feet
        body_length /= 3.281
        nose_length /= 3.281
        boattail_length /= 3.281
        body_diameter /= 3.281
        nose_diameter /= 3.281
        motor_diameter /= 3.281

        # change to handle approx zero machs too
        if mach == 0:
            return (0, 0, 0, 0)
        
        V = mach * speedOfSound
        q = (1.0/2.0) * density * V**2
        q_psf = q * 0.020885434273039 # Pascals to psf

        total_length = body_length + nose_length + boattail_length
        c_do_body_friction_missile = 0.053 * (total_length / body_diameter) * math.pow((mach/(q_psf * total_length)), 0.2)
        
        if (mach >= 1.0):
            c_do_body_wave_missile = (1.586 + 1.834/(mach**2)) * math.pow(math.atan(0.5/(nose_length/nose_diameter)), 1.69)
            if isPowered:
                c_do_body_base_missile = (1 - math.pow(motor_diameter, 2)/math.pow(body_diameter, 2)) * (0.25/mach)
            else:
                c_do_body_base_missile = 0.25/mach
        else:
            c_do_body_wave_missile = 0.0
            if isPowered:
                c_do_body_base_missile = (1 - math.pow(motor_diameter, 2)/math.pow(body_diameter, 2)) * (0.12 + 0.13 * mach**2)
            else:
                c_do_body_base_missile = 0.12 + 0.13 * mach**2

        c_do_body_missile = c_do_body_friction_missile + c_do_body_wave_missile + c_do_body_base_missile
        
        return (c_do_body_missile, c_do_body_friction_missile, c_do_body_wave_missile, c_do_body_base_missile)

    @staticmethod
    def c_do_fin(mach, speedOfSound, density, finArea, meanAreaChord, referenceArea, nFins):
        V = mach * speedOfSound
        q = (1./2.) * density * V**2
        q_psf = q * 0.020885434273039
        meanAreaChord_ft = meanAreaChord / 12.
        c_do_surface_friction = nFins * 0.0133 * math.pow(mach/(q_psf * meanAreaChord_ft), 0.2) \
                                * 2 * (finArea / referenceArea)
        
        return c_do_surface_friction

def coefficientOfFriction(Re, Re_crit = 500000):
    c_f_box = 0.
    
    if Re == 0:
        c_f_box = 0.
    elif Re < Re_crit:
        c_f_box = 1.328/math.sqrt(Re)
    else:
        beta = Re_crit * (0.074/math.pow(Re, 1.0/5.0) - 1.328/math.sqrt(Re))
        c_f_box = 0.074/math.pow(Re, 1.0/5.0) - beta/Re

    return c_f_box

class CoefficientsBox:
    def __init__(self):
        pass

    @staticmethod
    def c_do_body(mach, speedOfSound, density, kin_viscosity,\
        body_length, nose_length, boattail_length,\
        body_diameter, nose_diameter, boattail_diameter, Re_crit = 500000):

        if mach == 0:
            return (0, 0, 0, 0)

        V = mach * speedOfSound
        total_length = boattail_length + body_length + nose_length
        
        Re = V * total_length / kin_viscosity

        c_f_box = 0
        c_do_base_box = 0
        c_do_body_box = 0
        c_do_body_wave_missile = 0

        c_f_box = coefficientOfFriction(Re)

        c_do_body_box = (1 + 60.0/math.pow(total_length/body_diameter, 3.0) + 0.0025 * (body_length/body_diameter))\
        * (2.7*(nose_length/body_diameter) + 4*(body_length/body_diameter) + 2*(1-boattail_diameter/body_diameter)*(boattail_length/body_diameter))\
        * c_f_box

        c_do_base_box = 0.029 * math.pow(body_diameter/boattail_diameter, 3) / math.sqrt(c_do_body_box)

        # this is actually the correlation from Tactical Missile Design
        if (mach >= 1):
            c_do_body_wave_missile = (1.586 + 1.834/(mach**2)) * math.pow(math.atan(0.5/(nose_length/nose_diameter)), 1.69)
        else:
            c_do_body_wave_missile = 0

        c_do_body_box = compressibilityCorrection(c_do_body_box, mach)
        c_do_base_box = compressibilityCorrection(c_do_base_box, mach)
        c_do_total_box = c_do_body_box + c_do_base_box + c_do_body_wave_missile
        
        return (c_do_total_box, c_do_body_box, c_do_body_wave_missile, c_do_base_box)

    @staticmethod
    def c_do_fin_friction(mach, speedOfSound, density, kin_viscosity,\
        finThickness, nFins, bodyDiameter,\
        finRootChord, finTipChord, finSpan, finMidChord):

        if (mach == 0):
            return 0.

        V = mach * speedOfSound
        Re = density * V * finMidChord / kin_viscosity

        c_friction = coefficientOfFriction(Re)
        exposed_area = (1./2.) * (finRootChord + finTipChord) * finSpan
        planform_area  = exposed_area + (1./2.) * bodyDiameter * finRootChord

        c_do_fin_friction = 2 * c_friction * (1 + 2*(finThickness/finMidChord)) * (4*nFins*planform_area) / (math.pi * bodyDiameter**2)
        c_do_fin_friction = compressibilityCorrection(c_do_fin_friction, mach)
        return c_do_fin_friction

    # this function is from TMD
    # referenceArea is the body reference area!
    # angles should be in radians
    @staticmethod
    def c_do_fin_wave(mach, LESweepAngle, LETotalAngle, maxThickness,\
        finSpan, referenceArea, nFins):
        
        mach_le = mach * math.cos(LESweepAngle)
        # only good for 2 or 4 fins
        n_surface = nFins / 2
        
        if (mach_le < 1):
            return 0
        else:
            return ( n_surface * (1.429/mach_le**2)\
            * (math.pow(1.2*mach_le**2, 3.5) * 2.4/math.pow(2.8*mach_le**2 - 0.4, 2.5) - 1) \
            * math.pow(math.sin(LETotalAngle), 2) * math.cos(LESweepAngle) \
            * maxThickness * finSpan / referenceArea )

def getDragCoefficient(rocket, isPowered, atmosphere, mach, dragModel, constantDrag = 0):
    # the only key that is needed to be returned is 'all'
    # the others are optional
    results = {}
    
    if dragModel == DragModel.TacticleMissileDesign:
        # for now assume its always not powered
        c_d_all = CoefficientsTacticleMissile.c_do_body(
            mach, atmosphere.speedOfSound, atmosphere.density, isPowered,\
            rocket.geometry["body_length"], rocket.geometry["nose_length"], rocket.geometry["boattail_length"],\
            rocket.geometry["body_diameter"], rocket.geometry["nose_diameter"], rocket.geometry["motor_diameter"]
        )

        results['all'] = c_d_all[0]
        results['body_all'] = c_d_all[0]
        results['body_friction'] = c_d_all[1]
        results['body_wave'] = c_d_all[2]
        results['body_base'] = c_d_all[3]

    elif dragModel == DragModel.Constant:
        results['all'] = constantDrag
    elif dragModel == DragModel.CoefficientsBox:
        c_d_body = CoefficientsBox.c_do_body(
            mach, atmosphere.speedOfSound, atmosphere.density, atmosphere.kin_viscosity,\
            rocket.geometry["body_length"], rocket.geometry["nose_length"], rocket.geometry["boattail_length"],\
            rocket.geometry["body_diameter"], rocket.geometry["nose_diameter"], rocket.geometry["boattail_diameter"]
        )
        c_d_fin_friction = CoefficientsBox.c_do_fin_friction(
            mach, atmosphere.speedOfSound, atmosphere.density, atmosphere.kin_viscosity,\
            rocket.geometry["finThickness"], rocket.geometry["nFins"], rocket.geometry["body_diameter"],\
            rocket.geometry["finRootChord"], rocket.geometry["finTipChord"], rocket.geometry["finSpan"],\
            rocket.geometry["finMidChord"]
        )
        c_do_fin_wave = CoefficientsBox.c_do_fin_wave(
            mach, rocket.geometry["finLESweep"], rocket.geometry["finLETotal"],\
            rocket.geometry["finThickness"], rocket.geometry["finSpan"],\
            rocket.geometry["referenceArea"], rocket.geometry["nFins"]
        )

        c_d_all = c_d_body[0] + c_d_fin_friction + c_do_fin_wave
        c_d_fin_all = c_d_fin_friction + c_do_fin_wave

        results['all'] = c_d_all
        
        results['body_all'] = c_d_body[0]
        results['body_friction'] = c_d_body[1]
        results['body_wave'] = c_d_body[2]
        results['body_base'] = c_d_body[3]
        
        results['fin_all'] = c_d_fin_all
        results['fin_friction'] = c_d_fin_friction
        results['fin_wave'] = c_do_fin_wave

    return results