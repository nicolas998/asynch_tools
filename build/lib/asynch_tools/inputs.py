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
        #Read the topology 
        



