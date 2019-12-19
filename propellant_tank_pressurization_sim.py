''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * File name  :  RocketTrajectoryMain.py
 * Purpose    :  Assess propellant system required mass and tank volume
 * @author    :  Harrison Leece
 * Date       :  2019-12-20
 * Notes      :  None
 * Warnings   :  None
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Revision History
 * ================
 *   Ver      Date     Modified by:  Reason for change or modification
 *  -----  ----------  ------------  ---------------------------------------------------------------------
 *  1.0.0  2020-01-01  Harrison L    Initial release to LMARS
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '''

import numpy as np
import CoolProp.CoolProp as cp
import matplotlib.pyplot as plt

class GasSystem:
    '''
    This class encapsulated the propellant pressurization system used for the Loyola MARS
    rocket.  This includes the 4 important features:  Gas pressure tank, Gas Pressure regulator
    Fuel tank and Oxidizer tank.
    '''
    def __init__(self, prop_tank1, prop_tank2, supply_tank, step):
        self.step = step
        self.prop_tank1 = prop_tank1
        self.prop_tank2 = prop_tank2
        self.supply_tank = supply_tank

        '''
        Initialize a list for each tank which can store information used for
        plotting.
        '''
        #initialize lists used for plotting supply tank
        self.st_mass_list = []
        self.st_pressure_list = []
        self.st_temp_list = []
        #initialize lists used for plotting prop_tank1
        self.pt1_mass_list = []
        self.pt1_pressure_list = []
        self.pt1_temp_list = []
        self.pt1_propmass_list = []
        #initialize lists used for plotting prop_tank2
        self.pt2_mass_list = []
        self.pt2_pressure_list = []
        self.pt2_temp_list = []
        self.pt2_propmass_list = []
        #initialize some lists for plotting non-tank parameters
        self.time_list = []

    def simulate(self):
        t=0
        req_mass_flow = 0
        while ((self.supply_tank.p > self.prop_tank1.p_target) and (self.prop_tank1.prop_mass>0)):
            '''
            Append values to list BEFORE the simulation to capture inital values
            '''
            self.time_list.append(t)

            self.st_mass_list.append(self.supply_tank.mass)
            self.st_pressure_list.append(self.supply_tank.p)
            self.st_temp_list.append(self.supply_tank.temp)

            self.pt1_mass_list.append(self.prop_tank1.gas_mass)
            self.pt1_pressure_list.append(self.prop_tank1.p)
            self.pt1_temp_list.append(self.prop_tank1.temp)
            self.pt1_propmass_list.append(self.prop_tank1.prop_mass)

            self.pt2_mass_list.append(self.prop_tank2.gas_mass)
            self.pt2_pressure_list.append(self.prop_tank2.p)
            self.pt2_temp_list.append(self.prop_tank2.temp)
            self.pt2_propmass_list.append(self.prop_tank2.prop_mass)

            '''
            Actually runs the simulation
            '''
            t+=self.step
            #Calculates gas temperature out of tank using Joule Thomson effect
            exit_gas_temp = self.supply_tank.calc_regulator_exit_temp(self.prop_tank1.p_target)
            #Calculates mass flows into tank required to maintain the same specific volume
            #The specific volume at pressure and temperature is computed with CoolProp
            mass_in1 = self.prop_tank1.calc_req_mass(self.step, exit_gas_temp)
            mass_in2 = self.prop_tank2.calc_req_mass(self.step, exit_gas_temp)
            required_mass_flow = mass_in1 + mass_in2
            #Adiabatic expansion of the pressure supply tank
            self.supply_tank.tank_expansion(required_mass_flow, self.step)
            #add mass to the propellant tanks
            self.prop_tank1.update_ullage_mass_and_temp(mass_in1, exit_gas_temp)
            self.prop_tank2.update_ullage_mass_and_temp(mass_in2, exit_gas_temp)

    def plot_system(self):
        #supply tank plots
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.plot(self.time_list, self.st_pressure_list)
        plt.grid()
        plt.title('Pressure vs time')
        plt.ylabel('Pressure (Pa)')
        #velocity subplot
        plt.subplot(3, 1, 2)
        plt.title('Temperature vs Time')
        plt.ylabel('Temperature (Kelvin)')
        plt.plot(self.time_list, self.st_temp_list)
        plt.grid()
        plt.subplot(3, 1, 3)
        plt.title('Temperature vs Time')
        plt.ylabel('Mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.st_mass_list)
        plt.grid()

        #propellant tank 1 plots
        plt.figure()
        plt.subplot(4, 1, 1)
        plt.plot(self.time_list, self.pt1_pressure_list)
        plt.grid()
        plt.title('Pressure vs time')
        plt.ylabel('Pressure (Pa)')
        #velocity subplot
        plt.subplot(4, 1, 2)
        plt.title('Temperature vs Time')
        plt.ylabel('Temperature (Kelvin)')
        plt.plot(self.time_list, self.pt1_temp_list)
        plt.grid()
        plt.subplot(4, 1, 3)
        plt.title('Ullage Temperature vs Time')
        plt.ylabel('Mass (kg)')
        plt.plot(self.time_list, self.pt1_propmass_list)
        plt.grid()
        plt.subplot(4,1,4)
        plt.title('Prop Mass Remaining vs Time')
        plt.ylabel('Propellant mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.pt2_propmass_list)
        plt.grid()

        #propellant tank 2 plots
        plt.figure()
        plt.subplot(4, 1, 1)
        plt.plot(self.time_list, self.pt2_pressure_list)
        plt.grid()
        plt.title('Pressure vs time')
        plt.ylabel('Pressure (Pa)')
        #velocity subplot
        plt.subplot(4, 1, 2)
        plt.title('Ullage Temperature vs Time')
        plt.ylabel('Temperature (Kelvin)')
        plt.plot(self.time_list, self.pt2_temp_list)
        plt.grid()
        plt.subplot(4, 1, 3)
        plt.title('Temperature vs Time')
        plt.ylabel('Mass (kg)')
        plt.plot(self.time_list, self.pt2_mass_list)
        plt.grid()
        plt.subplot(4,1,4)
        plt.title('Prop Mass Remaining vs Time')
        plt.ylabel('Propellant Mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.pt2_propmass_list)
        plt.grid()

        plt.show()



class isochoricTank:
    '''
    This class describes the constant volume tank expelling intert gas to the
    two propellant tanks, regulated by a gas pressure regulator.
    The enthalpy drop over the gas pressure regulator is 0
    The Joule-Thomson effect is used to compute the exit temperature of the
    gas, which is then used to aid the assesment of how much mass is required
    out of the high pressure gas supply tank
    '''

    def __init__(self, p_i, t_i, volume, fluid_string):
        '''
        @parameters
        p_i: Initial pressure of the tank, in Pascals
        t_i: Initial temperature of the tank, in Kelvin
        volume: Initial volume of the tank.  This is a design quantity
        fluid_string: The cool-prop string which identifies the pure substance
                      in the tank.  "Helium" or "Nitrogen" typically
        '''
        #Set initial values for presssure and testing
        self.p_i = p_i
        self.temp_i = t_i
        #instance variables for pressure and temp
        self.p = p_i
        self.temp = t_i

        self.flu = fluid_string
        #initial density of the tank, kg/m^3
        self.rho_i = cp.PropsSI('D', 'P', p_i, 'T', t_i, self.flu)
        #instance variable for density
        self.rho = self.rho_i
        #volume in m^3.  This variable does not change with time
        self.volume = volume
        #initial mass in kg
        self.m_i = volume*self.rho_i
        #instance variable for mass
        self.mass = self.m_i

        self.update_specific_heats(p_i,t_i)
        #initialize lists used for plotting
        self.time_list = []
        self.mass_list = []
        self.pressure_list = []
        self.temp_list = []


    def update_specific_heats(self, pressure, temperature):
        '''
        Sets the instance variables for the specific heats and their ratio
        to the correct value for a given temperature and pressure with CoolProp
        '''
        self.c_p = cp.PropsSI('C', 'P', pressure, 'T', temperature, self.flu)
        self.c_v = cp.PropsSI('O', 'P', pressure, 'T', temperature, self.flu)
        self.gamma = self.c_p/self.c_v

    def calc_regulator_exit_temp(self, exit_pressure):
        '''
        Calculates the temperature of the gas leaving the  pressure regulator.
        This function captures the Joule-Thomson effect over the regulator.
        The Joule-Thomson effect will cause heating of Helium at high initial temperatures,
        but will cause cooling of Nitrogen and Helium in most cases.
        '''
        upstream_enthalpy = cp.PropsSI('H', 'P', self.p, 'T', self.temp, self.flu)
        downstream_temp = cp.PropsSI('T', 'H', upstream_enthalpy, 'P', exit_pressure, self.flu)
        return downstream_temp

    def tank_expansion(self, mass_out, step):
        '''
        Calculates the pressure and temperature of the gas in the tank as the
        gas expands.  This expansion is adiabatic and reversible.
        '''
        init_rho = self.rho_i
        #update mass
        self.mass = self.mass - mass_out
        new_rho = self.mass/self.volume

        init_pressure = self.p_i
        init_temperature = self.temp_i

        pressure = init_pressure*(new_rho/init_rho)**self.gamma
        temperature = init_temperature*(new_rho/init_rho)**(self.gamma - 1)
        #calculate new specific heats and gamma for the gas at the new
        #pressure and temperature
        self.update_specific_heats(self.p, self.temp)
        #Update the pressure and temperature of the system
        self.p = pressure
        self.temp = temperature

    def constant_mdot_test(self):
        '''
        Test function for the solving proceedure.
        Utilizes a constant mass flow rate out of the tank.
        Acts as a unit test for isochoric tanks, as the mass flow rate
        out of the tank is decoupled from the
        '''
        #step size in seconds.  Smaller numbers are more accurate
        step = .001
        #mass flow of gas out of tank.  Constan for testing
        mass = .00002
        t = 0
        while (self.p > 2930000):
            #Add initial states to the instance lists
            self.time_list.append(t)
            self.mass_list.append(self.mass)
            self.pressure_list.append(self.p)
            self.temp_list.append(self.temp)
            #time step
            t+=step
            #Find pressure and temperature in tank from ratio of density
            #mass/self.volume represents the new density of the tank after the
            #time step
            self.tank_expansion(mass, step)

    def plot_tank(self):
        plt.subplot(3, 1, 1)
        plt.plot(self.time_list, self.pressure_list)
        plt.grid()
        plt.title('Pressure vs time')
        plt.ylabel('Pressure (Pa)')
        #velocity subplot
        plt.subplot(3, 1, 2)
        plt.title('Temperature vs Time')
        plt.ylabel('Temperature (Kelvin)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.temp_list)
        plt.grid()
        plt.subplot(3, 1, 3)
        plt.title('Temperature vs Time')
        plt.ylabel('Mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.mass_list)
        plt.grid()
        plt.show()


class propellantTank:
    '''
    This object is for propellant tanks.  A target ullage pressure is known, and
    the ullage volume and gas temperature changes over the flight of the rocket.
    v_i is initial ullage pressure.  All inputs in SI
    '''
    def __init__(self, target_pressure, t_i, init_prop_mass, prop_mass_flow_rate, prop_density, initial_ullage_volume,fluid_string):
        self.flu = fluid_string
        self.rho_prop = prop_density
        self.ullage = initial_ullage_volume
        self.prop_volume = init_prop_mass/prop_density
        self.prop_mass = init_prop_mass
        #mass flow rate!
        self.prop_flow_rate = prop_mass_flow_rate
        #m^3/s of propellant required
        self.d_volume = prop_mass_flow_rate/prop_density
        self.p_target = target_pressure
        self.p = target_pressure
        self.temp = t_i
        self.gas_mass = initial_ullage_volume * cp.PropsSI('D', 'P', target_pressure, 'T', t_i, self.flu)
        #initialize lists used for plotting
        self.time_list = []
        self.prop_mass_list = []
        self.gas_mass_list = []
        self.gas_mass_input_list = []
        self.pressure_list = []
        self.temp_list = []

    def calculate_specific_heats(self, pressure, temperature):
        '''
        Sets the instance variables for the specific heats and their ratio
        to the correct value for a given temperature and pressure with CoolProp
        '''
        self.c_p = cp.PropsSI('C', 'P', pressure, 'T', temperature, self.flu)
        self.c_v = cp.PropsSI('O', 'P', pressure, 'T', temperature, self.flu)
        self.gamma = self.c_p/self.c_v

    def calc_req_mass(self, step, supply_temp):
        #Calcute volume change in tank per step
        volume_delta = self.d_volume * step
        #initialize variables outside of loop
        mass_required = 0
        m_guess = 0
        diff = 1
        while(np.abs(diff)>.001):
            m_guess += .0001
            guess_temp = self.compute_ullage_temp(m_guess, supply_temp)
            rho_guess = cp.PropsSI('D', 'P', self.p_target, 'T', guess_temp, self.flu)
            mass_required = volume_delta * rho_guess
            diff = mass_required - m_guess
        #update tank states
        self.gas_mass_input_list.append(mass_required)
        self.ullage = self.ullage + volume_delta
        self.prop_mass = self.prop_mass - self.prop_flow_rate * step
        return mass_required

    def update_ullage_mass_and_temp(self, added_gas_mass, added_gas_temp):
        self.gas_mass += added_gas_mass
        self.temp = self.compute_ullage_temp(added_gas_mass, added_gas_temp)

    def compute_ullage_temp(self, added_gas_mass, added_gas_temp):
        temp = (self.temp * self.gas_mass + added_gas_temp * added_gas_mass)/(self.gas_mass + added_gas_mass)
        return temp

    def proptank_test_no_mass(self):
        '''
        Acts as a unit test for isobaric tanks
        Utilizes a constant mass flow rate into the tank.
        '''
        #step size in seconds.  Smaller numbers are more accurate
        step = .005
        t = 0
        while (self.prop_mass > 0):
            #Add initial states to the instance lists
            self.time_list.append(t)
            #This should be the same number the whole time for this test
            self.gas_mass_list.append(self.gas_mass)
            self.prop_mass_list.append(self.prop_mass)
            self.pressure_list.append(self.p)
            self.temp_list.append(self.temp)
            #time step
            t+=step
            self.prop_mass = self.prop_mass - self.prop_flow_rate * step

            self.temp = temp
            self.p = p

    def plot_tank(self):
        plt.subplot(4, 1, 1)
        plt.plot(self.time_list, self.pressure_list)
        plt.grid()
        plt.title('Pressure vs time')
        plt.ylabel('Pressure (Pa)')
        #velocity subplot
        plt.subplot(4, 1, 2)
        plt.title('Temperature vs Time')
        plt.ylabel('Temperature (Kelvin)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.temp_list)
        plt.grid()
        plt.subplot(4, 1, 3)
        plt.title('Temperature vs Time')
        plt.ylabel('Propellant mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.prop_mass_list)
        plt.grid()
        plt.subplot(4, 1, 4)
        plt.title('Gas input vs Time')
        plt.ylabel('Gas mass (kg)')
        plt.xlabel('Time (s)')
        plt.plot(self.time_list, self.gas_mass_input_list)
        plt.grid()
        plt.show()


#If running the program from this file, the below code executes
if __name__ == "__main__":
    o_f_ratio = 1.9
    total_prop_flow_rate = 5.739
    fuel_flow_rate = total_prop_flow_rate/(o_f_ratio+1)
    ox_flow_rate = total_prop_flow_rate-fuel_flow_rate

    rho = cp.PropsSI('D', 'T', 298.15, 'P', 101325, 'Nitrogen')
    print("Density of Nitrogen at STP: ", rho)

    he_tank = isochoricTank(34470000,298,.03,"Helium")
    print("Helium initial mass", he_tank.m_i, 'kg')

    lox_tank = propellantTank(3103000, 298, 127.79, ox_flow_rate, 1146, .003, "Helium")
    jeta_tank = propellantTank(3103000, 298, 67.25, fuel_flow_rate, 819, .003, "Helium")

    #he_tank.constant_mdot_test()
    #he_tank.plot_tank()

    #print('Lox tank initial propellant mass: ', 127.79, 'kg')
    #print('Lox tank propellant mass flow rate: ', lox_tank.prop_flow_rate,'kg/s')

    gas_system = GasSystem(lox_tank, jeta_tank, he_tank, .001)
    gas_system.simulate()
    gas_system.plot_system()
