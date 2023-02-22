#!/usr/bin/env python

"""
L2Pred
LambdaPred

"""

import copy, time, warnings
import scipy as sp
import numpy as np
from scipy import linalg
from collections import OrderedDict, defaultdict


class BasePred():
    
    def _checktype(self, obj, classnames):
        if not (type(obj).__name__ in list(classnames)): raise TypeError(f'{type(obj)} not allowed as linkdata input. Must be {classnames}')        
    
    def _order(self, input_lst):
        assert type(input_lst) is list
        if self.shuffle:
            return input_lst[np.random.permutation(len(input_lst))]
        else:
            return input_lst
        
    def _return_stddev(self, _stansda):
        _stddev = _stansda.stats[:,1][:,np.newaxis]
        if np.any(np.isinf(_stddev)): warnings.warn('inf value detected in standardizer. implies maf=0')
        _stddev[np.isinf(_stddev)] = 1
        return _stddev
    
    def predict(self, *, srd, n_inchunk=1000, stansda=None):
        if self.verbose: print('predicting for given srd:')
        beta = np.array(self.beta_est).flatten()
        if stansda is None:
            stansda = self.stansda
        if (stansda is None) and (self.stansda is None)  :
            stansda = self.linkdata.get_combined_unit_stansda()
        assert len(beta) == srd.shape[1]
        start=0; y_pred_lst = [];
        for j in range(0, srd.shape[1]+n_inchunk, n_inchunk)[1:]:
            stop = min(j, srd.shape[1])
            print('chunk of start:stop=',start, stop, end='\r') if self.verbose else None
            slice_srd = srd[:,start:stop]
            sda = slice_srd.read().standardize(stansda)
            assert np.all(sda.sid == slice_srd.sid)
            y_pred_lst.append(sda.val.dot(beta[start:stop]))
            start = stop
        y_pred = np.vstack(y_pred_lst).sum(axis=0)
        self.iid_pred = copy.copy(srd.iid)
        self.y_pred = y_pred
        return y_pred
    

class L2Pred(BasePred):
    
    def __init__(self, linkdata, *, h2, n_iter=1, 
                 random_state=None, shuffle=False,
                 local_rm=False, 
                 clear_linkdata=True,
                 verbose=False):
        
        # Stuff all the args into fields.
        _excl_lst = ['self', 'kwg_dt']
        kwg_dt = {key: item for key, item in locals().items() if not (key in _excl_lst)}
        for key, item in locals().items():
            if not (key in _excl_lst): 
                self.__setattr__(key, item)

        # Checks:
        self._checktype(linkdata, 'LinkageData')
        assert type(h2) is float
    
    def get_weights(self, return_frame=True):
        if not return_frame: raise NotImplementedError('contact dev') 
        weights = self._sst_df.copy()
        weights['weight'] = self._weights
        return weights
    
    def _compute_cur_beta_tilde(self, *, beta, i_reg, linkdata):
        
        beta_tilde = linkdata.get_beta_marginal_region(i=i_reg)
        
        if self.local_rm:
            # RM
            True
            
        return beta_tilde
    
    
    def fit(self): return self.sample()
    
    def sample(self):
        
        s=self; 
        tpl=s.n_iter, s.verbose, s.h2, s.linkdata # Could have a code reading function here..
        (     n_iter,   verbose,   h2,   linkdata) = tpl
        
        beta_mrg = linkdata.get_beta_marginal()
        p        = len(beta_mrg)
        sst_df   = linkdata.get_sumstats_cur()
        n_eff    = sst_df['n_eff'].median()
        i_lst    = linkdata.get_i_list()
        assert linkdata.n_snps_total == p
        
        # Random state & Shuffling:
        if self.random_state != None:
            np.random.seed(self.random_state)
    
        # Initalisations:
        beta  = np.zeros((p,1))
        
        # Sampling Loops:
        if verbose: print('Starting iterations of Sampler:')
        for itr in range(n_iter):
            
            for i_reg in self._order(i_lst):
                do_show = ((i_reg % 10 == 0) | (i_reg<10)) & verbose
                if do_show: print(f'-> itr={itr}, i_reg={i_reg} <-  ', end='\r')

                # Old logistical stuff:
                D = linkdata.get_auto_linkage_region(i=i_reg)
                p_reg = len(D); assert p_reg > 0 # I think this is not needed anymore.
                idx_reg = range(*linkdata.get_auto_range_region(i=i_reg))
                
                # RM:
                beta_tilde = self._compute_cur_beta_tilde(beta=beta, 
                                                                i_reg=i_reg, linkdata=linkdata)

                # Sample beta from MVN:
                tau = p*(1-h2)/(n_eff*h2)
                beta_i = linalg.solve(D + tau*np.eye(p_reg), beta_tilde)
                beta[idx_reg] = beta_i
                
        
        # Post proc & storage:
        self.stansda = linkdata.get_stansda()
        _stddev = self._return_stddev(self.stansda)
        self._weights = beta/_stddev
        self._stddev  = _stddev
        self._sst_df  = sst_df
        if self.clear_linkdata: del self.linkdata

        if verbose: print('----- Done with Sampling -----')
            
        return self

class LambdaPred(BasePred):
    pass
