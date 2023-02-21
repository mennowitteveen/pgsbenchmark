#!/usr/bin/env python

"""
PrivacyPreservingMetricsComputer
MultiPGSComputer

"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from pysnptools.standardizer import UnitTrained
# from functools import partial
try:
    from tqdm.auto import tqdm
except:
    tqdm = lambda x:x
    
# tqdm = partial(tqdm, position=0, leave=True)
# from tqdm import tqdm
# from functools import partial
# tqdm = partial(tqdm, position=0, leave=True)
# tqdm = partial(tqdm, position=0, leave=True, ncols=70, delay=0.3)

locals_dt = dict()    


class PrivacyPreservingMetricsComputer():
    
    def __init__(self, *, linkdata, brd, Bm, s=None, dtype='float32',
                 clear_linkage=True, pbar=tqdm, verbose=True):
        
        self.linkdata   = linkdata
        self.brd        = brd
#         assert (np.isnan(s).sum()+np.isinf(s).sum()) == 0
        self.s          = s
        self.Bm         = Bm 
        self.dtype      = dtype
        self.clear_linkage = clear_linkage
        if not pbar: pbar = lambda x: x
        self.pbar       = pbar
        self.verbose    = verbose
            
    def evaluate(self, debug=False):
        
        # Load and init variables:
        linkdata = self.linkdata
        #linkdata = self.linkdata.init() # init, in case required.
        brd = self.brd; Bm = self.Bm
        if self.verbose & (self.s is None): print('Retrieving Standard Dev. var \'s\'')
        s = linkdata.get_stansda().stats[:,[1]] if self.s is None else self.s
        assert (np.isnan(s).sum()+np.isinf(s).sum()) == 0 # add s =standard dev as argument with object creation if this line keeps failing
        bCb = 0.; BmBt = 0.
        info_dt = dict()

        # Cycle through the blocks:
        for i, geno_dt in self.pbar(linkdata.reg_dt.items()):
            #if self.verbose: print(f'PPB: Processing region {i}', end='\r')
            # Ready the LD:
            L = linkdata.get_left_linkage_region(i=i)
            D = linkdata.get_auto_linkage_region(i=i)
            R = linkdata.get_right_linkage_region(i=i)
            lr = linkdata.get_left_range_region(i=i)
            ar = linkdata.get_auto_range_region(i=i)
            rr = linkdata.get_right_range_region(i=i)

            # Ready The Weights:
            B_L = brd[:,lr[0]:lr[1]].read().val.astype(self.dtype).T
            B_D = brd[:,ar[0]:ar[1]].read(dtype=self.dtype).val.T
            B_R = brd[:,rr[0]:rr[1]].read(dtype=self.dtype).val.T
            B_L = s[lr[0]:lr[1]]*B_L
            B_D = s[ar[0]:ar[1]]*B_D
            B_R = s[rr[0]:rr[1]]*B_R

            # Do the computation:
            CB = L.dot(B_L) + D.dot(B_D) + R.dot(B_R)
            bCb += (B_D*CB).sum(axis=0)
            BmBt += (B_D.T.dot(Bm.iloc[ar[0]:ar[1],:])).T
            info_dt[i] = dict(shapeL=L.shape, shapeD=D.shape, shapeR=R.shape, 
                              lr=lr, ar=ar, rr=rr)

            # Pruning to minimize memory overhead:
            if (i > 0) and self.clear_linkage:
                linkdata.clear_linkage_region(i=i-1)
            #if (i > 38) & debug: break; return locals()

        # Complete resutls:
        linkdata.clear_all_xda()
        cols = brd.row.astype(str).flatten()
        bCb  = pd.DataFrame(bCb[np.newaxis,:], index=['bCb'], columns=cols)
        BmBt = pd.DataFrame(BmBt, index=Bm.columns, columns=cols)
        ppbr2_df = (BmBt**2)/bCb.loc['bCb']
        res_dt = dict(ppbr2_df=ppbr2_df, bCb=bCb, BmBt=BmBt, info_dt=info_dt, s=s)

        return res_dt
    
# locals_dt = dict()
class MultiPGSComputer():
    
    def __init__(self, *, brd, unscaled=True, verbose=False, dtype='float32', allow_nan=False, pbar=tqdm):
        self.brd   = brd
        if hasattr(brd, 'val'):
            assert np.sum(np.isnan(brd.val)) == 0 
        self.unscaled = unscaled
        self.verbose = verbose
        self.dtype  = dtype
        self.allow_nan = allow_nan
        self.pbar = pbar
        
    def predict(self, *, srd, prd=None, n_inchunk=1000, stansda=None):
            
        # Load that PGS (& optionaly phenos)
        brd = self.brd
        Yhat = np.zeros((srd.shape[0], brd.shape[0]), dtype=self.dtype)        
        assert np.all(brd.col.astype(str) == srd.sid)
        if prd: 
            assert np.all(srd.iid == prd.iid)
            pda   = prd.read(dtype=self.dtype).standardize()
            Ytru  = pda.val
            Bm    = np.zeros((srd.shape[1], Ytru.shape[1])) + np.nan
        
        # Loop through Genome:
        stansda_lst = []; start=0
        for start in self.pbar(range(0, srd.shape[1], n_inchunk)):
            stop = min(start+n_inchunk, srd.shape[1])
            sda, stansda = srd[:,start:stop].read(dtype=self.dtype).standardize(return_trained=True)
            X = sda.val
            s = stansda.stats[:,1][:,np.newaxis]; s[np.isinf(s)] = 1
            L = brd[:,start:stop].read(dtype=self.dtype).val
            B = s*L.T # s seems of little effect on time here, projected loading takes time 200s for HM3 8K betas. for 10k induv.
            Yhat += X@B
            if prd: Bm[start:stop] = X.T@Ytru
            stansda_lst.append(stansda)
            
        if prd:    
            Bm = Bm/Ytru.shape[0]
            if not self.allow_nan: assert np.isnan(Bm).sum() == 0
            if not self.allow_nan: assert np.isnan(Bm).sum() == 0
            Bm = pd.DataFrame(Bm, index=srd.sid, columns=prd.col)
            Ytru = pd.DataFrame(Ytru, # Make Ytru a proper dataframe
                index=pd.MultiIndex.from_arrays(prd.iid.T, names=('fid','iid')),
                columns=prd.col)
        else:
            Ytru=None; Bm=None
            
        # Combine Standardizers:
        sid     = np.concatenate([stan.sid   for stan in stansda_lst])
        assert  np.unique(sid).shape[0] == sid.shape[0]
        stats   = np.concatenate([stan.stats for stan in stansda_lst])
        stansda = UnitTrained(sid, stats)   
        s = stansda.stats[:,1][:,np.newaxis]; s[np.isinf(s)] = 1
        
        # Create Yhat dataframe:
        Yhat  = pd.DataFrame(
            data    = Yhat, 
            index   = pd.MultiIndex.from_arrays(srd.iid.T, names=('fid','iid')),
            columns = self.brd.row.astype(str)
        ); assert Yhat.isna().sum().sum() == 0
        
        res_dt = dict(Yhat=Yhat, Bm=Bm, brd=brd, Ytru=Ytru, stansda=stansda, s=s)
        
        return locals()
    
    def run(self):
        pass
    
    def fit(self):
        pass
