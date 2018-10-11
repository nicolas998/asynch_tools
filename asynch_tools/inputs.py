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


#Basin topo class used to generate input rainfall and initial conditions for
#asynch
class basinTopo:

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
        Area = []; Slope=[]; Long = []
        for l in L[3::3]:
            Param = [float(i) for i in l.split()]
            Area.append(Param[0])
            Slope.append(Param[1])
            Long.append(Param[2])
        return 'vacio'

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
                __read_basin_prm_file__()
            else:
                print('Waning: not a valid .prm file extension')
        #Read lookup file 
        if path_lookup is not None:
            name, ext = os.path.splitext(path_lookup)
            if ext == '.lookup':
                self.path_lookup = path_lookup
                __read_basin_lookup_file__()
            else:
                print('Warning: not a valid .lookup file')

    def __RainUniformRainfall__(N, u_min, u_max):
        '''Generate uniform random rainfall'''
        return np.random.uniform(u_min, u_max, N)

    def __RainHistogramRainfall__(Histogram, Ngen = 50):
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

    def __RainWriteVariableHill__(f,rain, timeStep, Hill, maskHills):
        '''Updates variable hill rainfall file'''
        f.write('%d\n' % Hill)
        f.write('%d\n' % rain.size)
        for c,r in enumerate(rain):
            Step = (c+1)*timeStep
            f.write('%.3f %.3f \n' % (Step, r))
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
            MaskedHills = HillsIDs < HillsIDs.max()+1

        #Llama funciones para generar la lluvia
        if VariableHills:
            #Opens file and writes the total number of hills in it.
            f = open(path, 'w')
            f.write('%d \n\n\n' % HillsIDs.size)
            #Writes rainfall for each hill
            for Hill, Mask in zip(HillsIDs, MaskedHills):
                #Writes if hill is not masked
                if Mask:
                    #Selects rainfall function 
                    if rain == 'urandom':
                        # Uniform rainfall function
                        rainValues = __RainUniformRainfall__(Nrain, urand_min, urand_max)
                        #Writes rainfall for hill
                        __RainWriteVariableHill__(f,rainValues, RainTimeStep, Hill, MaskedHills)
                    elif rain == 'hrandom':
                        # Histogram based rainfall 
                        rainValues = __RainHistogramRainfall__(BaseHietogram, Ngen)
                        #Writes rainfall for hill
                        __RainWriteVariableHill__(f,rainValues, RainTimeStep, Hill, MaskedHills)
            #Close rainfall file 
            f.close()
        elif VariableHills is False:
            #Opens file and writes the total number of time steps.
            f = open(path, 'w')
            f.write('%d \n\n' % Nrain)
            #Selects rainfall depending on the generator function
            if rain == 'urandom':
                rainValues = __RainUniformRainfall__(Nrain, urand_min, urand_max)
            #Writes rainfall for each time step 
            for c,R in enumerate(rainValues):
                Step = (c+1)*RainTimeStep
                f.write('%.3f %.3f \n' % (Step, R))
            #Close rainfall file 
            f.close()









       




