# Something python bash python something

# Doing something to make all the import suggestions in the code etc work nicely will come later
# this requires writing fancy __init__.py files and having subdirs, _stuff.py files and more..

# I put this here, because I had a hunch it would generate desired python behavior.
# __all__ = ['Tree','fullvars','sizegb','beep','Timer','implot', 'redo_all_above', 'do_all_above', 'Struct'] 

import time
import os
import ntpath
import pickle
import inspect
from sys import getsizeof
from collections import defaultdict
from IPython.display import Audio, display, Javascript
from IPython import get_ipython
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd

# Tree:
def Tree():
    return defaultdict(Tree)

#fullvars:
def fullvars(obj):
    return {key: getattr(obj, key) for key in dir(obj)}

# Size computer:
def sizegb(var, verbose=False):
    size_in_gb = getsizeof(var)/1024**3
    if verbose:
        print('GB: ', size_in_gb)
    return size_in_gb

# Beeper:
def beep():
    framerate = 4410
    play_time_seconds = 1
    audio_data = [np.sin(2*np.pi*300*t)+np.cos(2*np.pi*240*t) 
                  for t in np.linspace(0, play_time_seconds, framerate*play_time_seconds)]
    return Audio(audio_data, rate=framerate, autoplay=True)

# Timer:
class Timer():
    def __init__(self):
        self.time_dt = {}
        self.timeepoch_cnt = -1

    def tic(self, dashline=True):
        if dashline:
            print('-'*42)
        self.timeepoch_cnt += 1
        cnt = self.timeepoch_cnt
        self.time_dt[cnt] = [time.time()]
        self.time_lst = self.time_dt[cnt]

    def toc(self, rep=''):
        cnt = self.timeepoch_cnt
        time_lst = self.time_lst
        t = time.time()
        current_line = inspect.currentframe().f_back.f_lineno
        print(f'E{cnt} Step({len(time_lst)}): {t-time_lst[0]:.3f}, {t-time_lst[-1]:.3f}  <-- line({current_line}) -- {rep}')
        time_lst.append(t)
        

def find_varname(var):
    lcls = inspect.stack()[2][0].f_locals
    for name in lcls:
        if id(var) == id(lcls[name]):
            return name
    return None

# Dumping this here, burr hurr
# def find_varname_better(a, f, b): # needs to be inside..
#     frame = inspect.currentframe()
#     frame = inspect.getouterframes(frame)[1]
#     string = inspect.getframeinfo(frame[0]).code_context[0].strip()
#     args = string[string.find('(') + 1:-1].split(',')
    
#     names = []
#     for i in args:
#         if i.find('=') != -1:
#             names.append(i.split('=')[1].strip())

#         else:
#             names.append(i)
    
#     print(names)

# import re
# import traceback

# def func(var, otherarg):
#     stack = traceback.extract_stack()
#     filename, lineno, function_name, code = stack[-2]
#     vars_name = re.compile(r'\((.*?)\).*$').search(code).groups()[0]
#     print(vars_name)
#     return

# foobar = 5
# durr = 3
# func(foobar + 5, 6)
        
# Implot:
def implot(M, cmap=None, title=None):
    if title is None:
        #title = find_varname(M) # so so..
        frame = inspect.currentframe()
        frame = inspect.getouterframes(frame)[1]
        string = inspect.getframeinfo(frame[0]).code_context[0].strip()
        args = string[string.find('(') + 1:-1].split(',')
#         if args[-1][-1] == ')': args[-1] = args[-1][:-1]
        names = []
        for i in args:
            if i.find('=') != -1:
                names.append(i.split('=')[1].strip())
            else:
                names.append(i)    
        title = names[0]
        
    plt.imshow(M, cmap=cmap, interpolation='none'); plt.colorbar(); plt.title(title); plt.show()
    
# Redo all above:
def redo_all_above():
    display(Javascript('Jupyter.notebook.kernel.restart(); IPython.notebook.execute_cells_above();'))
    # Best 2 run in middle of browser page, wait and visual view will end up in roughly same place.
    
# Redo all above:
def do_all_above():
    display(Javascript('IPython.notebook.execute_cells_above();'))
    # Best 2 run in middle of browser page, wait and visual view will end up in roughly same place.
    
class Struct(dict):
    def __init__(self, dt=None):
        """Convert a dictionary to a class
        @param :adict Dictionary
        """
        if dt is not None:
            self.__dict__.update(dt)
                
    def update(self, dt):
        """Convert a dictionary to a class
        @param :adict Dictionary
        """
        assert type(dt) is dict
        super().update(dt)
        self.__dict__.update(dt)
        
def psrc(obj, return_source=False):
    """Print the code of a Python object
    """
    src = inspect.getsource(obj)
    print(src)
    if return_source:
        return src

def jobinfo(job, return_string=False):
    string = get_ipython().system('jobinfo {job.job_id}')
    print(string)
    if return_string:
        return string
    
    
def corr(X, Y=None, ddof=0): # so a correlation might not actually have a ddof.
    if Y is None: Y = X 
    assert type(X) is pd.DataFrame
    assert X.shape[0] == Y.shape[0]
    X = (X-X.mean())/X.std(ddof=ddof)
    Y = (Y-Y.mean())/Y.std(ddof=ddof)
    arr = X.values.T.dot(Y.values)/(X.shape[0]-ddof)
    C = pd.DataFrame(arr, index=X.columns, columns=Y.columns)
    return C

def pcorr(X, Y, ignore_nans=False):

    if np.isnan(X).any()  or np.isnan(Y).any():
        if not ignore_nans:
            raise Exception('Input contains NaNs, this is not allowed')
            
    n = X.shape[0]

    m_X = np.mean(X, axis=0)[np.newaxis,:]
    m_Y = np.mean(Y, axis=0)[np.newaxis,:]
    # sn_X = np.sum( ,axis=0)
    # sn_Y = np.sum( ,axis=0)
    s_X = np.std(X, axis=0)[np.newaxis,:]
    s_Y = np.std(Y, axis=0)[np.newaxis,:]
    s_XY = s_X.T.dot(s_Y)

    R = (X-m_X).T.dot(Y-m_Y)
    R = (R/s_XY)/n

    # Compute P-values matrix:
    dist = sp.stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
    P = 2*dist.cdf(-abs(R))
    
    return R, P

class ResultsStorage():

    def __init__(self, *, timer, suffix='xp', res_base_dn='./results/', verbose=True):

        self.timer=timer
        self.timer.tic(dashline=False) # Start timer.
        self.suffix = suffix
        self.res_base_dn = res_base_dn
        self.ts = pd.Timestamp.now()
        self.ts_str = self.ts.floor('s').isoformat().replace('T','_')
        self.res_dn = os.path.join(res_base_dn, f'{self.ts_str}_{suffix}/')
        self.res_dn = os.path.abspath(os.path.expanduser(self.res_dn))
        if verbose:
            print('ResultStorage dir : ', self.res_dn)
        assert not os.path.exists(self.res_dn)
        os.makedirs(self.res_dn)


    def write_pickle(self, sub_dn='', mode='wb', **kwargs):
        assert len(kwargs) == 1
        key = list(kwargs.keys())[0]
        fn = os.path.join(self.res_dn, sub_dn, key + '.pickle')
        fn = self._prepare_path(fn)
        with open(fn, mode) as f:
            pickle.dump(kwargs[key], f)

    def write_nbformat(self, nb, *, fn, mode='w'):
        fn = self._prepare_path(os.path.join(self.res_dn, fn))
        with open(os.path.join(self.res_dn, fn), mode) as f:
            nbformat.write(nb, f)

    def write_string(self, string, *, fn, mode='w'):
        fn = self._prepare_path(os.path.join(self.res_dn, fn))
        with open(fn, mode) as f:
            f.write(string)

    def _prepare_path(self, fn):
        full_dn = ntpath.dirname(fn)
        os.makedirs(full_dn, exist_ok=True)
        return fn

    def report_completion(self):

        t = self.timer.time_dt[0][0]
        dt = str(time.time()-t) + 's'
        report_string = f"""
        Completed!
        Total-time-taken: {str(pd.Timedelta(dt))}
        """
        self.write_string(report_string, fn='completed.txt')

        