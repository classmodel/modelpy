# class4gl extends the standard 'model' environment to be able to take global air profiles as input 

from model import model as class4gl
from model import model_output as class4gl_output
import numpy as np

class air_input(object):
    def __init__(self):
        self.status = 'init'
        
        
    def set_air_input(self,INPUT):
        PARAMS,ONE_COLUMN = INPUT.PARAMS,INPUT.ONE_COLUMN

        self.PARAMS = PARAMS
        self.ONE_COLUMN = ONE_COLUMN
        self.__dict__.update(PARAMS.to_dict()['value'])
        self.__dict__.update(ONE_COLUMN.to_dict('list'))
        self.status = 'filled'
        
        #convert all list to arrays for CLASS
        for key in self.__dict__.keys():
            if type(self.__dict__[key]).__name__ == 'list':
                self.__dict__[key] = np.array(self.__dict__[key])

from model import model_input
class4gl_input = type('class4gl_input', (model_input,air_input), dict(c='c'))

