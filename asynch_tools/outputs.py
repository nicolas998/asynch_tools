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


def __ReadStagesNamesfromGlobal__(path):
    '''Reads the names and order of the stages defined in the global file'''
    #Check extension
    name, ext = os.path.splitext(path) 
    if ext != '.gbl':
        print('Error: not a .gbl file')
        return 1
    #opens the file and read lines
    f = open(path, 'r')
    Lines = f.readlines()
    f.close()
    #Read the names
    pos = Lines.index('%Components to print\n')
    NStages = int(Lines[pos+1])
    Names = []
    for i in range(pos+2,pos+2+NStages):
        Names.append(Lines[i].split()[0])
    return Names

def ReadDatFile(path, gbl_path = None):
    '''Given the path it reads the .dat file of an asynch simulation'''
    #Check extension
    name, ext = os.path.splitext(path) 
    if ext != '.dat':
        print('Error: not a .dat file')
        return 1
    #opens the file and read lines
    f = open(path,'r')
    Lines = f.readlines()
    f.close()
    #Read file parameters
    Ncols = int(Lines[1])
    NrecHills = int(Lines[0])
    NumRec = int(Lines[3].split(' ')[1])
    HillsSim = [int(i.split(' ')[0]) for i in Lines[3::NumRec+2]]
    #Read data
    Coef = 0
    Suma = 4
    Data = []
    for i in range(NrecHills):
        D = []
        for l in Lines[Suma+NumRec*Coef:NumRec*(Coef+1)+Suma]:
            D.append([float(i) for i in l.split(' ')[:-1]])
        Data.append(np.array(D).T)
        Coef += 1
        Suma += 2
    #Transform data into a pandas DataFrame
    Col1 = []
    Col2 = []
    for h in HillsSim:
        for i in range(Ncols-1):
            Col1.append(str(h))
            Col2.append('Stage'+str(i))
    Columns = [Col1, Col2]
    B = np.array([i[1:] for i in Data])
    index = Data[0][0]
    A = B.reshape(B.shape[0]*B.shape[1], B.shape[2])
    return pd.DataFrame(A.T, index=index, columns=Columns)