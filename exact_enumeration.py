import matplotlib.pyplot as plt
import numpy as np 
from math import sqrt 
from math import floor
from math import ceil
from math import exp 
import pickle
from operator import itemgetter
from matplotlib import cm, colors
from numpy import amin, amax, ravel
from scipy.optimize import curve_fit
from scipy.special import logsumexp


def get_all_saws(current_paths, length, current_conformation):
    """ Все созданные конформации лежат в первом аргументе функции get_all_saws, 
    поэтому первым аргументом нужно подавать переменную, которой ранее присвоили пустой список
    
    """
    if(length==1):
        current_paths.append(current_conformation)
    else:
        for step in [(1, 0), (-1, 0), (0, 1),  (0, -1)]:
            new_point = (current_conformation[-1][0]+step[0], current_conformation[-1][1]+step[1] )
            if new_point in current_conformation:
                continue
            else:
                temp_path = current_conformation.copy()
                temp_path.append(new_point)
                get_all_saws(current_paths, length-1,  temp_path)
                
                
                
def get_sequences(length):
    """Получаем все последовательности заданной длины, состоящие из 0 и 1.
    1 - это H
    0 - это Р"""
    if(length ==1):
        return [[0], [1]]
    else:
        previous = get_sequences(length - 1)
        result = [] 
        for i in range(len(previous)):
            current = previous[i].copy()
            current.append(0)
            result.append(current.copy())
            current[-1] = 1
            result.append(current)
            #print(result) 
        return sorted(result)
        
 
        
class Protein(object):
    def __init__(self, sequence, conformation):
        self.sequence = sequence
        self.conformation = conformation
        
    def count_proteins_contacts(self):
        """Для данного белка считает число топологических контактов HH и контакты HP/PP.
        Возвращает кортеж из двух элементов: первый элемент - число контактов НН, второй элемент - остальные контакты."""
        hh=0
        hp_pp = 0
        steps = [(1, 0), (-1, 0), (0, 1),  (0, -1)]
        for i in range(1, len(self.conformation)-1):
            not_topological=[self.conformation [i-1], self.conformation[i+1]]
            for step in steps:
                new_point = (self.conformation[i][0]+step[0], self.conformation [i][1]+step[1] )
                if (new_point in self.conformation  and (new_point not in not_topological)):
                    position = self.conformation.index(new_point)
                    if( self.sequence[position]==1 and self.sequence[i]==1 ):
                        hh=hh+1
                    else:
                        hp_pp=hp_pp + 1 
        for step in steps:
            new_point_begin = (self.conformation[0][0]+step[0], self.conformation[0][1]+step[1] )
            new_point_end= ( self.conformation[-1][0]+step[0], self.conformation [-1][1]+step[1])
            if(new_point_begin in self.conformation  and new_point_begin!=self.conformation [1]):
                position = self.conformation.index(new_point_begin)
                if( self.sequence[position]==1 and self.sequence[0]==1 ):
                    hh=hh+1
                else:
                    hp_pp=hp_pp + 1 
            if( new_point_end in self.conformation  and new_point_end!=self.conformation[-2]):
                position = self.conformation.index(new_point_end) 
                if( self.sequence[position]==1 and self.sequence[-1]==1 ):
                    hh=hh+1
                else:
                    hp_pp=hp_pp + 1 
        return (hh//2, hp_pp//2)
    
    def count_x(self): 
        steps = [(1, 0), (-1, 0), (0, 1),  (0, -1)]
        x=[]
        for i in range(len(self.conformation)):
            k=0
            for step in steps: 
                point = ( self.conformation[i][0]+step[0], self.conformation [i][1]+step[1]   )
                if(point in self.conformation):
                    k=k+1
                else:
                    break
            if(k==4):
   
                x.append(self.sequence[i])
                
        return x 
        
 

