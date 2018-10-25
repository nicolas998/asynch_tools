#!Tools to generate input information for asynch runs
#!Copyright (C) <2016>  <Nicolas Velasquez Giron>

#!This program is free software: you can redistribute it and/or modify
#!it under the terms of the GNU General Public License as published by
#!the Free Software Foundation, either version 3 of the License, or
#!(at your option) any later version.

#!This program is distributed in the hope that it will be useful,
#!but WITHOUT ANY WARRANTY; without even the implied warranty of
#!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#!GNU General Public License for more details.

#!You should have received a copy of the GNU General Public License
#!along with this program.  If not, see <http://www.gnu.org/licenses/>.

#External packages
import numpy as np
import os
import pandas as pd

#GlobalGenerator: class used to generate global files for asynch, with this 
#you can specify output files, time deltas, running time, etc.
class GlobalGenerator:
    
    def __UpdateGlobalVariables__(self, key):
        '''Function to update the global variables in the file'''
        #Variables with similar behavior 
        ListSimilar = ['Topology', 'Parameters','Initial','Rainfall','Links2Save']
        #Temporal dictionary and position in the global file 
        Dtemp = self.GlobalVariables[key]
        p = self.editGlobal.index(Dtemp['Phrase'])
        #Finds position in the global file 
        if key == 'Time':
            for i,param in zip(Dtemp['Pos'], Dtemp['param']):
                self.editGlobal[p+i] = '%.1f\n' % param
        elif key == 'Globals':
            Dtemp = self.GlobalVariables[key]
            Texto = '%d' % len(Dtemp['param'])
            pos = Dtemp['Pos']
            for i in Dtemp['param']:
                if np.abs(i)>0.01 or i == 0.0:
                    Texto += ' %.4f' % i
                else:
                    Texto += ' %.3E' % i
            self.editGlobal[p+pos] = Texto + '\n'
        elif key == 'Out_hydro' and self.GlobalVariables[key]['param'] is not None:
            Dtemp = self.GlobalVariables[key]
            pos = Dtemp['Pos']
            Texto = str(Dtemp['TimeDelta']) + ' ' + Dtemp['param'] + '\n'
            self.editGlobal[p+pos] = self.editGlobal[p+pos][0] + ' ' + Texto
        else:
            try:
                KeyInPos = ListSimilar.index(key)
                if self.GlobalVariables[key]['param'] is not None:
                    Dtemp = self.GlobalVariables[key]
                    pos = Dtemp['Pos']
                    self.editGlobal[p+pos] = self.editGlobal[p+pos][0] + ' ' + Dtemp['param']+'\n'
            except:
                pass
    
    def __init__(self, baseGlobalPath = None):
        '''Global generator for asynch runs.'''
        self.GlobalVariables = {'Time':{'param': [0, 20*24*3600],
                'Phrase': '%Begin and end date time\n',
                'Pos': [1, 2]},
            'Globals':{'param': [0.33, 0.20, -0.1, 0.0, 2.0425e-6, 0.02, 0.5, 0.10, 0.0, 99.0, 3.0, 0.75],
                'Phrase': '%Global parameters\n',
                'Pos': 2},
            'Topology':{'param': None,
                'Phrase': '%Topology (0 = .rvr, 1 = database)\n',
                'Pos': 1},
            'Parameters': {'param': None,
                'Phrase': '%DEM Parameters (0 = .prm, 1 = database)\n',
                'Pos': 1},
            'Initial': {'param': None,
                'Phrase': '%Initial state (0 = .ini, 1 = .uini, 2 = .rec, 3 = .dbc)\n',
                'Pos': 1},
            'Rainfall': {'param': None,
                'Phrase': '%Rain\n',
                'Pos': 1},
            'Links2Save': {'param': None,
                'Phrase': '%.sav files for hydrographs and peak file\n',
                'Pos': 2},
            'Out_hydro': {'param': None,
                'Phrase': '%Where to put write hydrographs\n',
                'Pos': 2,
                'TimeDelta': 10},
                }
        if baseGlobalPath is not None:
            #Opens the base global and puts it as a back variable
            f = open(baseGlobalPath, 'r')
            self.baseGlobal = f.readlines()
            f.close()
            #Make the editable global
            self.editGlobal = self.baseGlobal.copy() 

    def ChangeVariable(self, VarName, VarValue):
        '''Function to change the value of a variable for the global file to write'''
        self.GlobalVariables[VarName]['param'] = VarValue
        self.__UpdateGlobalVariables__(VarName)
        
    def WriteGlobalFile(self, GlobalPathName, Return2BaseGlobal = False):
        '''Function to write a global file with the changes done'''
        #Writes new global file into disk
        f = open(GlobalPathName, 'w')
        f.writelines(self.editGlobal)
        f.close()
        #If yes, return the editGlobal to baseGlobal
        if Return2BaseGlobal:
            self.editGlobal = self.baseGlobal.copy()
    
    
#Basin topo class used to generate input rainfall and initial conditions for
#asynch
class basinTopo:

    ########################################################################################################
    #Functions read structures and to initialize basin object 
    
    def __read_basin_from_msg__(self):
        '''Read hillsIds and hillsConect from a msg file'''
        #Read the data frame in pandas format and obtains the FullStruct
        self.FullStruct = pd.read_msgpack(self.path)
        #Obtains the hillsIds
        self.hillsIds = np.array([int(i) for i in self.FullStruct.index.tolist()])
        #Obtains the hillsConnect
        self.hillsConnect = self.FullStruct['Parents'].to_dict()

    def __read_basin_from_rvr__(self):
        '''Read hillsIds and HillsConnect from a .rvr file
        do not produce a FullStruct file'''
        #Read the .rvr file 
        f = open(self.path,'r')
        L = f.readlines()
        f.close()
        #Obtains the hills IDs
        self.hillsIds = np.array([int(i) for i in L[2::3]])
        #Obtaines the hillsConnect
        self.hillsConnect = {}
        for ids, par in zip(L[2::3],L[3::3]):
            Parents = par.split()
            ids = ids.split()[0]
            if len(Parents) == 1:
                self.hillsConnect.update({ids: ()})
            else:
                self.hillsConnect.update({ids: tuple([int(i) for i in par.split()[1:]])})

    def __read_basin_prm_file__(self):
        '''Read a parameter file from a .prm'''
        #Read
        f = open(self.path_prm,'r')
        L = f.readlines()
        f.close()
        #Obtain parameters
        ParamsDict = {'p1':[],'p2':[],'p3':[]}
        for l in L[3::3]:
            Param = [float(i) for i in l.split()]
            for i in range(3):
                ParamsDict['p'+str(i+1)].append(Param[i])
        for i in range(3):
            ParamsDict['p'+str(i+1)] = np.array(ParamsDict['p'+str(i+1)])
        self.params = ParamsDict

    def __read_basin_lookup_file__(self):
        return 'vacio'

    def __init__(self, path, path_prm = None, path_lookup = None):
        '''basinTopo: element with the topological structure
        of a basin that can be run by asynch.

        Parameters:
            - path: path to the .rvr file containing the topology
                of the basin
            - path_prm: path to the .prm file with the properties.
            - path_lookup: path to the .lookup file with the coordinates'''
        #Set the path of the base file and read hills and connect
        self.path = path
        name, ext = os.path.splitext(path)
        if ext == '.msg':
            self.__read_basin_from_msg__()
        elif ext == '.rvr':
            self.__read_basin_from_rvr__()
        #Read parameters file
        if path_prm is not None:
            name, ext = os.path.splitext(path_prm)
            if ext == '.prm':
                self.path_prm = path_prm
                self.__read_basin_prm_file__()
            else:
                print('Waning: not a valid .prm file extension')
        #Read lookup file 
        if path_lookup is not None:
            name, ext = os.path.splitext(path_lookup)
            if ext == '.lookup':
                self.path_lookup = path_lookup
                self.__read_basin_lookup_file__()
            else:
                print('Warning: not a valid .lookup file')

    ########################################################################################################
    #Functions to obtain random rainfal fields 
                
    def __RainUniformRainfall__(self, N, u_min, u_max):
        '''Generate uniform random rainfall'''
        return np.random.uniform(u_min, u_max, N)

    def __RainHistogramRainfall__(self, Histogram, Ngen = 50):
        '''Generate rainfall copying a base histogram corresponding
        to mean intensity values during a storm
        Parameters:
            - Histogram: Base histogram (hietogram)
            - Ngen: Number of generations used to copy the hietogram
                large values imply very similar rainfall, low numbers 
                imply coarse rainfall'''
        #Make the histogram a float so can use to generate numbers
        Histogram = Histogram.astype(float)
        Histogram = np.append(0.0, Histogram)
        #Estimate probability, and copy it 
        Prob = Histogram / Histogram.sum()
        hist,bins = np.histogram(np.random.uniform(0,1,Ngen), bins=Prob.cumsum())
        hist = hist.astype(float) / hist.sum()
        #Returns rainfall estimated
        return hist * Histogram.sum()

    def __RainWriteVariableHill__(self, f,rain, timeStep, Hill, maskHills):
        '''Updates variable hill rainfall file'''
        f.write('%d \n' % Hill)
        f.write('%s\n' % str(rain.size+2))
        f.write('0 0.00\n')
        for c,r in enumerate(rain):
            Step = (c+1)*timeStep
            f.write('%.3f %.3f \n' % (Step, r))
        f.write('%d 0.0\n' % ( (c+2)*timeStep) )
        f.write('\n')
    
    def Rain2strFile(self, path, rain = 'urandom', Nrain = None, RainTimeStep = 5,
        VariableHills = False, urand_min = 0, urand_max = 10, MaskedHills = None,
        BaseHietogram = None, Ngen = 50):
        '''Code to obtain an .str file with the shape:

        Structure for variable Hills file (VariableHills = True):

        {number of links}
        {link id 1} {number of changes}
        {time 1} {value 1}
        {time 2} {value 2}
        {time 3} {value 3}

        Structure for uniform Hills file (VariableHills = False):

        {number of changes}
        {time 1} {value 1}
        {time 2} {value 2}
        {time 3} {value 3}

        more info see:
        https://asynch.readthedocs.io/en/latest/input_output.html#storm-files

        Parameters:
            - path: Path to save rainfall .str
            - HillsIDs: Vector with the IDs of the hills.
            - rain: Rainfall generator function
                urandom: Uniform random generator (urand_min, urand_max)
                hrandom: Histogram based random generator
                    - BaseHietogram: Base hietogram used to obtain random rainfall similar to it
                    - Ngen: Number of generations used to imitate BaseHietogram (50).
                        High values: this imply that generated hietograms are similar to BaseHietogram
                        Low values: Imply a coarser hietogram (higher intensities) 
                np.array: copia un array de Nhills x Nrain
            - RainTimeStep: Time step between rainfall intervals (minutes)
            - Nrain: Number of rainfall inputs.
            - VariableHills (False): Obtains same or different rainfall for each
                hill
            - MaskedHills: Boolean vector which hills has rainfall'''

        #Type of extension
        extension = '.ustr'
        if VariableHills: extension = '.str'
        #Path Fix
        name, ext = os.path.splitext(path) 
        if ext != extension:
            path = name + extension

        #Mascara de los hills 
        if MaskedHills is None:
            MaskedHills = self.hillsIds < self.hillsIds.max()+1

        #Llama funciones para generar la lluvia
        if VariableHills:
            #Opens file and writes the total number of hills in it.
            f = open(path, 'w')
            f.write('%d \n\n\n' % self.hillsIds.size)
            #Writes rainfall for each hill
            for Hill, Mask in zip(self.hillsIds, MaskedHills):
                #Writes if hill is not masked
                if Mask:
                    #Selects rainfall function 
                    if rain == 'urandom':
                        # Uniform rainfall function
                        rainValues = self.__RainUniformRainfall__(Nrain, urand_min, urand_max)
                        #Writes rainfall for hill
                        self.__RainWriteVariableHill__(f,rainValues, RainTimeStep, Hill, MaskedHills)
                    elif rain == 'hrandom':
                        # Histogram based rainfall 
                        rainValues = self.__RainHistogramRainfall__(BaseHietogram, Ngen)
                        #Writes rainfall for hill
                        self.__RainWriteVariableHill__(f,rainValues, RainTimeStep, Hill, MaskedHills)
                else:
                    #Rain values equal to zero 
                    if rain == 'urandom': rainValues = np.zeros(Nrain)
                    if rain == 'hrandom': rainValues = np.zeros(BaseHietogram.size)
                    self.__RainWriteVariableHill__(f,rainValues, RainTimeStep, Hill, MaskedHills)
            #Close rainfall file 
            f.close()
            return rainValues
        elif VariableHills is False:
            #Selects rainfall depending on the generator function
            if rain == 'urandom':
                #Uniform random generator
                rainValues = self.__RainUniformRainfall__(Nrain, urand_min, urand_max)
            if rain == 'hrandom':
                #Hietogram random generator
                rainValues = self.__RainHistogramRainfall__(BaseHietogram, Ngen)
                Nrain = rainValues.size
            #Opens file and writes the total number of time steps.
            f = open(path, 'w')
            f.write('%d \n\n' % Nrain)
            #Writes rainfall for each time step 
            for c,R in enumerate(rainValues):
                Step = (c+1)*RainTimeStep
                f.write('%.3f %.3f \n' % (Step, R))
            #Close rainfall file 
            f.close()
            return rainValues
        
    ########################################################################################################
    #Functions to obtain random initial conditions 
            
    def __InitialStatesUrandom__(self, rmin, rmax):
        '''Generate a initial state uniform random vector rmax and rmin
        correspond to the min and max values of the random generator'''
        #Generator 
        States = []
        for rmi,rma in zip(rmin, rmax):
            States.append(np.random.uniform(rmi,rma))
        return States

    def __InitialStatesNrandom__(self, mean, desv):
        '''Generate a initial state normal random vector mean and dev
        correspond to the parameters of the normal distribution,
        the function is truncated in order to obtain possitive values only'''
        #Generator 
        States = []
        for m,s in zip(mean, desv):
            Val = -9
            while Val<0:
                Val = np.random.normal(m,s)
            States.append(Val)
        return States

    def __InitialStatesWrite__(self, file, HillsID, states):
        '''Function to write ini files for asynch, requires 
        a file handler object, the hillID, and the vector with the states'''
        # Writes number of the link
        file.write('%d \n' % HillsID)
        for state in states:
            file.write('%.5f ' % state)
        file.write('\n')

    def Initial2IniFile(self, path, modelType = 252, initial = 'urandom', InitialTime = 0,
        VariableHills = False, rand_param1 = None, rand_param2 = None, **kwargs):
        '''Writes a file with the initial states of the hills, this file could be 
        .uini or .ini.

        ini: Different initial values for each hill: 

        {model type}
        {number of links}
        {initial time}
        {link id 1}
        {initial value 1} {initial value 2}
        {link id 2}
        {initial value 1} {initial value 2}

        uini: Same initial values for all hills:

        {model type}
        {initial time}
        {initial value 1} {initial value 2}

        more info see:
        https://asynch.readthedocs.io/en/latest/input_output.html#initial-values-input

        Parameters:
            - path: Path to save the .ini or .uini file.
            - modelType: Configuration of the model (default 252 Toplayer model)
            - HillsIDs: Number od the id for each hill apply only for variable hills (.ini)
            - initial: Method to obtain initial values.
                - urandom: Put a uniform random value at each link
                    rand_param1: minimum value of the uniform distribution
                    rand_param2: maximum value of the uniform distribution
                - nrandom: Put a normal random value at each link
                    rand_param1: mean value of the normal distribution
                    rand_param2: deviation value of the normal distribution
                - np.array: Put a variable based on an array containing stages.
            - Ntages: Number of stages for the model.  For information about this see:
                https://asynch.readthedocs.io/en/latest/builtin_models.html
            - InitialTime: Place to put the forcing into the simulation.
            - urand_min: List with the minimum initial values for random generation.
            - urand_max: List with the maximum initial values for random generation.
        '''

        #Type of extension
        if VariableHills:
            extension = '.ini'
        else:
            extension = '.uini'
        #Path Fix
        name, ext = os.path.splitext(path) 
        if ext != extension:
            path = name + extension

        #Void function just to compare if initial argument is a function
        VoidF = lambda x: 1

        if VariableHills:
            #Write a .ini file with different initial states for each hill
            f = open(path, 'w')
            f.write('%d \n' % modelType)
            f.write('%d \n' % self.hillsIds.size)
            f.write('%.5f \n' % InitialTime)
            #Iterate in al hills in order to write its initial values
            for Hill in self.hillsIds:
                #Selec initial state function 
                if initial == 'urandom':
                    #Estado inicial aleatorio uniforme 
                    States = self.__InitialStatesUrandom__(rand_param1, rand_param2)
                elif initial == 'nrandom':
                    #Estado inicial aleatorio uniforme 
                    States = self.__InitialStatesNrandom__(rand_param1, rand_param2)
                elif type(initial) == type(VoidF):
                    #Estado inicial dado por una funcion
                    States = initial(rand_param1, rand_param2)
                elif initial == np.ndarray:
                    #Estado inicial dado por un array con los valores
                    States = np.copy(initial[0])
                #Writes states into the ini file    
                self.__InitialStatesWrite__(f, Hill, States)
            #Close file 
            f.close()
        else:
            #Selec initial state function 
            if initial == 'urandom':
                #Estado inicial aleatorio uniforme 
                States = self.__InitialStatesUrandom__(rand_param1, rand_param2)
            elif initial == 'nrandom':
                #Estado inicial aleatorio uniforme 
                States = self.__InitialStatesNrandom__(rand_param1, rand_param2)
            elif type(initial) == type(VoidF):
                #Estado inicial dado por una funcion
                States = initial(rand_param1, rand_param2)
            elif initial == np.ndarray:
                #Estado inicial dado por un array con los valores
                States = np.copy(initial)
            #Writes a .uini file with the same states for all the hills 
            f = open(path, 'w')
            f.write('%d \n' % modelType)
            f.write('%.5f \n\n' % InitialTime)
            for state in States:
                f.write('%.5f ' % state)
            f.write('\n')
            f.close()








       




