#!/usr/bin/env python

"""
LinkageData
durr tst
"""
import scipy as sp
import numpy as np
import pandas as pd
from scipy import linalg
from sys import getsizeof

import warnings, importlib, json, os, glob
from collections import OrderedDict, deque, defaultdict
from pysnptools.standardizer import Unit, UnitTrained
import pysnptools as pst
# import pysnptools.util as pstutil
# from pysnptools.standardizer import UnitTrained
# from dataclasses import dataclass

class SqrtNinv(Unit):
    def __init__(self):
        super(SqrtNinv, self).__init__()

class BaseLinkageData():

    def __init__(self, *, sst_df=None, regdef_df=None, master_dt=None, #There should be sst_df or master_dt
                 srd=None, sda_standardizer=Unit,
                 prd=None, pda_standardizer=Unit,
                 lrd=None, lda_standardizer=None,
                 grd=None, gda_standardizer=False,
                 
                 shift=0, cm=None, _setzero=True, ddof=0, #ddof should remain 0 for now
                 
                 clear_xda=True, # Refactor with _clear?
                 clear_linkage=False,
                 compute_sumstats=False,
                 calc_allelefreq=False,
                 intersect_apply=True,
                 
                 _onthefly_retrieval=True, # These underscore options are the advanced developer options
                 _save_vars = ['L','D','R','sst_df'],
                 _clear_vars = ['L','D','R'],
                 _cross_chrom_ld = False,
                 _save_s2sst = True,
                 
                 gb_size_limit=10., dtype='float32', verbose=False):
        
        if True:
            # bim and fam df have to be supplied because pysnptools halvely
            # implemented these portions of the genetic data into their object
            # meaning that srd cannot be relied uppon
            excl_lst = ['self','kwg_dt','excl_lst']
            kwg_dt = {key: item for key, item in locals().items() if not (key in excl_lst)}
            for key, item in locals().items():
                if not (key in excl_lst): 
                    self.__setattr__(key, item)
            # New rule: blx have to be created from the inside
            # Perhaps later it can be made into a special load instead of a compute

            # first-checks & inits:
            if cm is not None: assert cm > 0
            if lrd is not None: raise NotImplementedError('lrd not possible atm.')
            if grd is not None:
                assert gda_standardizer or (gda_standardizer is None)
            assert type(compute_sumstats) is bool
            if ddof != 0: raise NotImplementedError('delta deg. of freedom have to 0 for this version')
            self.reg_dt = dict()
            self.cur_total_size_in_gb = 0.0
            self.xda_q = deque()
            [self.xda_q.append((-1,'')) for _ in range(5)]  # put 5x -1 in queue
            self.reloaded_xda_cnt = 0
            self._fn_lst = []

            # Checks            
            if srd is not None:
                assert type(sst_df) is pd.DataFrame
                self._check_xrd()
                assert isinstance(sst_df, pd.DataFrame)
                assert isinstance(regdef_df, pd.DataFrame)
                self.init_regions()
            elif master_dt is not None:
                # Fill attributes in case master_dt is present:
                for key, item in master_dt.items():
                    setattr(self, key, item)
                reg_dt=dict()
                for pre_i, geno_dt in self.reg_dt.items(): reg_dt[int(pre_i)] = geno_dt
                self.reg_dt = reg_dt # An ugly type conversion hack cause json does not allow i to be integer, but forces it to be a string.
            else:
                raise Exception('Essentials not present')

    def _check_xrd(self):

        if self.srd is not None:
            assert pst.snpreader.SnpReader in self.srd.__class__.__mro__

        if self.prd is not None:
            n_start = len(self.prd.iid)
            n_pheno = self.prd.shape[1]
            if n_pheno > 1: raise NotImplementedError(f'only one pheno in prd allowed {n_pheno} detected ({prd.col}), for now. remove other phenos')
            self.srd, self.prd = pst.util.intersect_apply([self.srd, self.prd])
            if len(self.prd.iid) != n_start:
                warnings.warn('Number of samples do not match up after internal intersection, samples were lost:' 
                              f'{n_start - len(self.prd.iid)}, start = {n_start}, after_intersection = {len(self.prd.iid)}')
                if not self.intersect_apply: raise Exception('Intersection was required, but may not performed. Hence raising this error.')

        if self.grd is not None:
            # Check alignment for now, auto alignment needs work cause iid stuffs:
            if self.srd is not None:
                if not np.all(self.grd.sid == self.srd.sid):
                    raise Exception('snps of grd and srd not matching up, align first,'
                                    ' auto align will be implemented later')
            else:
                raise NotImplementedError('Not sure what to do with grd if no srd is present. not implemented.')
        
    ###########################
    # Regions Administration:
    if True:

        def init_regions(self):
            do_beta_moving = ('beta_mrg' in self.sst_df.columns)
            if not do_beta_moving:
                warnings.warn('No \'beta_mrg\' column detected in sst_df! This means that no summary stats were detected.')
            else:
                assert 'n_eff' in self.sst_df.columns
            cur_chrom = None
            i = 0; n_snps_cumsum = 0
            sst_df_lst = []
            for reg_cnt, (_, row) in enumerate(self.regdef_df.iterrows()):
                # Move region into specialized dictionary
                regid = row['regid'];
                chrom = row['chrom']
                start = row['start'];
                stop  = row['stop']

                # Map Variants to region
                ind = self.sst_df.chrom == chrom
                ind = (self.sst_df['pos'] >= start) & ind
                ind = (self.sst_df['pos'] < stop) & ind
                sid = self.sst_df['snp'][ind].values
                indices = self.srd.sid_to_index(sid)  # if sid not strickly present this will give an error!
                n_snps_reg = len(indices)
                if n_snps_reg == 0:
                    continue
                else:
                    geno_dt = dict(regid=regid,
                                   chrom=chrom,
                                   start=start,
                                   stop=stop,
                                   start_j=n_snps_cumsum)
                    n_snps_cumsum += n_snps_reg
                    geno_dt['stop_j'] = n_snps_cumsum
                    sst_df = self.sst_df[ind].copy(); sst_df['i'] = i
                    geno_dt['sst_df'] = sst_df
                    assert geno_dt['start_j'] == sst_df.index[0]; sst_df_lst.append(sst_df)
                    assert geno_dt['stop_j']  == sst_df.index[-1] + 1
                    if do_beta_moving:
                        geno_dt['beta_mrg'] = geno_dt['sst_df']['beta_mrg'].copy().values[:, np.newaxis]
                        assert len(geno_dt['beta_mrg'].shape) == 2
                    if self.srd is not None:
                        geno_dt['srd'] = self.srd[:, indices]
                        geno_dt['stansda'] = self.sda_standardizer() if self.sda_standardizer is not None else None
                    else:
                        raise NotImplementedError()
                    if self.grd is not None:
                        geno_dt['grd'] = self.grd[:, indices]
                        geno_dt['stangda'] = self.gda_standardizer() if self.gda_standardizer is not None else None
                    # Count up if things are actually stored in reg_dt
                    self.reg_dt[i] = geno_dt
                    i += 1
            self.n_snps_total = n_snps_cumsum
            sst_df = pd.concat(sst_df_lst, axis=0)
            self.sst_df = sst_df

        def get_i_list(self):
            return list(self.reg_dt.keys())

        def _load_all_snpdata(self):
            # load all regions
            for i, geno_dt in self.reg_dt.items():
                sda = geno_dt['srd'].read(dtype=self.dtype)
                stansda = sda.train_standardizer(apply_in_place=True,
                                                 standardizer=geno_dt['stansda'])
                geno_dt['sda'] = sda
                geno_dt['stansda'] = stansda

    ###########################
    ## Compute: ###############

    # Local Linkage Stuff: ####
    if True:
    
        def compute_linkage_sameregion(self, *, i):
            return self.compute_linkage_shiftregion(i=i, shift=0)

        def regions_compatible(self, *, i, j):
            try:
                if self.reg_dt[i]['chrom'] == self.reg_dt[j]['chrom']:
                    res = True
                elif self._cross_chrom_ld:
                    res = True
                else:
                    res = False
            except Exception as e:
                if (not (i in self.reg_dt.keys())) or (not (j in self.reg_dt.keys())):
                    res = False
                else:
                    raise e
            return res

        def compute_linkage_shiftregion(self, *, i, shift):
            j = i + shift
            if self.regions_compatible(i=i, j=j):
                self_sda = self.get_sda(i=i)
                dist_sda = self.get_sda(i=j)
                n = len(self_sda.iid)
                S_shift = self_sda.val.T.dot(dist_sda.val) / (n - self.ddof)
                return S_shift
            else:
                self_sda = self.get_sda(i=i)
                return np.zeros((self_sda.val.shape[1], 0))

        def compute_linkage_cmfromregion(self, *, i, cm):
            geno_dt = self.reg_dt[i]; lst = []
            if cm < 0: # Doing left:
                stop_j   = geno_dt['start_j']
                cm_left  = geno_dt['sst_df'].loc[stop_j]['cm'] 
                slc_df = self.sst_df.loc[:stop_j-1]
                slc_df = slc_df[slc_df.chrom==geno_dt['chrom']]
                slc_df = slc_df[slc_df.cm > (cm_left + cm)]
                start_i = slc_df['i'].min()
                start_i = -7 if np.isnan(start_i) else start_i
                for cur_i in range(start_i, i):
                    lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_i-i))
                    if start_i == -7: break
                L = np.concatenate(lst, axis=1)[:,-slc_df.shape[0]:] # concat & clip
                if self._setzero:
                    cms_reg    = geno_dt['sst_df']['cm'].values
                    cms_distal = slc_df['cm'].values
                    cms_L      =  cms_distal[np.newaxis,:] - cms_reg[:,np.newaxis]
                    setzero_L  = cms_L < cm
                    L[setzero_L] = 0
                    assert L.shape == setzero_L.shape
                return L
            else:
                start_j   = geno_dt['stop_j']
                cm_right  = geno_dt['sst_df'].loc[start_j-1]['cm']
                slc_df = self.sst_df.loc[start_j:]
                slc_df = slc_df[slc_df.chrom==geno_dt['chrom']]
                slc_df = slc_df[slc_df.cm < (cm_right + cm)]
                stop_i = slc_df['i'].max()
                stop_i = i+2 if np.isnan(stop_i) else stop_i + 1
                for cur_i in range(i+1, stop_i):
                    lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_i-i))
                R = np.concatenate(lst, axis=1)[:,:slc_df.shape[0]] # concat & clip
                if self._setzero:
                    cms_reg    = geno_dt['sst_df']['cm'].values
                    cms_distal = slc_df['cm'].values
                    cms_R     =  cms_distal[np.newaxis,:] - cms_reg[:,np.newaxis]
                    setzero_R = cms_R > cm
                    R[setzero_R] = 0
                    assert R.shape == setzero_R.shape
                return R
        
    # Misc Stuff: #############
    if True:
    
        def compute_sumstats_region(self, *, i):
            geno_dt = self.reg_dt[i]
            sda = self.get_sda(i=i)
            X = sda.val
            y = self.get_pda().val
            n = len(y)
            c_reg = X.T.dot(y) / (n - self.ddof)
            return c_reg   

        def compute_allelefreq_region(self, *, i):
            # Speed might be improved by using dot prod here, instead of sums
            # np.unique was way slower (5x)
            geno_dt = self.reg_dt[i]
            n, p_blk = sda.val.shape
            sst_df = geno_dt['sst_df'].copy()
            cnt0   = np.sum(sda.val==0, axis=0)
            cnt1   = np.sum(sda.val==1, axis=0)
            cnt2   = np.sum(sda.val==2, axis=0)
            cntnan = np.sum(np.isnan(sda.val), axis=0)
            assert np.allclose(cnt0 + cnt1 + cnt2 + cntnan, n)
            sst_df['altcnt=0']   = cnt0
            sst_df['altcnt=1']   = cnt1
            sst_df['altcnt=2']   = cnt2
            sst_df['altcnt=nan'] = cntnan
            sst_df['altfreq']    = (cnt1 + cnt2)/(n - cntnan)
            sst_df['missfreq']   = 1 - cntnan/n
            return sst_df

        def compute_ldscores_region(self, *, i):
            sst_df = self.reg_dt[i]['sst_df'].copy()
            L = self.get_left_linkage_region(i=i)
            D = self.get_auto_linkage_region(i=i)
            R = self.get_right_linkage_region(i=i)
            for k, j in enumerate(sst_df.index):
                slds = np.sum(L[k]**2) + np.sum(D[k]**2) + np.sum(R[k]**2)
                sst_df.loc[j, 'lds'] = np.sqrt(slds)
            return sst_df
        
    ############################
    ## Retrieve: ###############
    
    # Local Linkage: ############
    if True:
    
        def retrieve_linkage_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_linkage_region(i=i)
            if self.verbose:   print('\nDone')
            if self.clear_xda: self.clear_all_xda()

        def retrieve_linkage_region(self, *, i):

            geno_dt = self.reg_dt[i]
            if 'store_dt' in geno_dt.keys():
                self.load_linkage_region(i=i)
            shift = self.shift; cm = self.cm
            compute_sumstats = self.compute_sumstats

            if 'L' in geno_dt.keys():
                if 'D' in geno_dt.keys():
                    if 'R' in geno_dt.keys():
                        return None  # everything is done now.

            if self.verbose: print(f'Computing LD for region #{i} on chr{geno_dt["chrom"]}', end='\r')
            # Refactor: if linkage is only in blocks this code will lead to recomputation...
            if (shift > 0):
                L_lst = []
                R_lst = []
                for cur_shift in range(1, shift + 1):
                    L_lst.append(self.compute_linkage_shiftregion(i=i, shift=-cur_shift))
                    R_lst.append(self.compute_linkage_shiftregion(i=i, shift=cur_shift))

                # Store Linkage in geno_dt
                geno_dt['L'] = np.concatenate(L_lst[::-1], axis=1)  # L stands for left
                geno_dt['D'] = self.compute_linkage_sameregion(i=i)  # Linkage within region, D is convention from LDpred 1
                geno_dt['R'] = np.concatenate(R_lst, axis=1)  # R stands for right

                # Indices needed for slicing and dicing matched variables (e.g. beta weights):
                geno_dt['start_j_L'] = geno_dt['start_j'] - geno_dt['L'].shape[1]
                geno_dt['stop_j_L'] = geno_dt['start_j']
                geno_dt['start_j_R'] = geno_dt['stop_j']
                geno_dt['stop_j_R'] = geno_dt['stop_j'] + geno_dt['R'].shape[1]

            elif (shift==0) and (cm is None):  # Only same region has to be done.
                geno_dt['D'] = self.compute_linkage_sameregion(i=i)

            elif (shift==0) and cm > 0:
                geno_dt['L'] = self.compute_linkage_cmfromregion(i=i, cm=-cm)
                geno_dt['D'] = self.compute_linkage_sameregion(i=i)
                geno_dt['R'] = self.compute_linkage_cmfromregion(i=i, cm=cm)

                # Indices needed for slicing and dicing matched variables (e.g. beta weights):
                geno_dt['start_j_L'] = geno_dt['start_j'] - geno_dt['L'].shape[1]
                geno_dt['stop_j_L'] = geno_dt['start_j']
                geno_dt['start_j_R'] = geno_dt['stop_j']
                geno_dt['stop_j_R'] = geno_dt['stop_j'] + geno_dt['R'].shape[1]

            if compute_sumstats:
                self.retrieve_sumstats_region(i=i)
              
        def load_linkage_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.load_linkage_region(i=i)
            if self.verbose: print('\nDone')
            
        def load_linkage_region(self, *, i):
            geno_dt = self.reg_dt[i]
            store_dt = geno_dt['store_dt']
            for varname, file_dt in store_dt.items():
                module = importlib.import_module('.'.join(file_dt['typestr'].split('.')[:-1]))
                cname  = file_dt['typestr'].split('.')[-1]
                CurClass = getattr(module, cname) # Retrieves module.submodule.submodule.. etc
                curfullfn = os.path.join(self.curdn, file_dt['fn'])
                geno_dt[varname] = CurClass(pd.read_hdf(curfullfn, key=file_dt['key']))
                if self.verbose: print(f'loading: fn={curfullfn} key={file_dt["key"]}'+' '*50, end='\r')
                
        def save(self, fn, keyfmt='ld/chrom{chrom}/i{i}/{varname}', fmt='hdf5', mkdir=False, dn=None):
            self.curdn = os.path.dirname(fn) if (dn is None) else dn
            fn = os.path.basename(fn) if (dn is None) else fn
            if mkdir: os.makedirs(self.curdn, exist_ok=True)
            if (fmt != 'hdf5'): raise Exception(f'Only hdf5 file format supported atm, not {fmt}') 
            for i, geno_dt in self.reg_dt.items():
                self.save_linkage_region(i=i, fn=fn)
                
            # Saving of 'logistical' data for the object
            master_lst = [ 'shift', 'cm', '_setzero',
             'clear_xda', 'clear_linkage', 'compute_sumstats', 'calc_allelefreq', 
             '_onthefly_retrieval', '_save_vars', '_clear_vars', 
             'gb_size_limit', 'dtype', 'verbose', 'n_snps_total']
            geno_lst = ['regid','chrom','start','stop','start_j','stop_j',
                        'start_j_L', 'stop_j_L', 'start_j_R', 'stop_j_R','store_dt']
                
            def caster(arg, types):
                if np.isscalar(arg):
                    if isinstance(arg, np.integer): arg = int(arg)
                if type(arg) is int: return int(arg)
                assert type(arg) in types
                return arg

#             if hasattr(self,'s'): assert (self.s.shape == (self.n_snps_total,1))
            master_dt = dict(); maxlen = 20
            for key in master_lst:
                var = getattr(self, key)
                if type(var) is list:
                    for item in var:
                        assert type(item) in (bool, str, float, int, type(None))
                        if type(item) is str: assert (len(item) < maxlen)
                elif type(var) is str:
                        assert len(var) < maxlen
                master_dt[key] = caster(var, (list, bool, float, int, str, type(None)))

            reg_dt = dict()
            for i, geno_dt in self.reg_dt.items():
                newgeno_dt = dict()
                for key in geno_lst:
                    if not (key in geno_dt.keys()): continue
                    newgeno_dt[key] = caster(geno_dt[key], (str, int, dict))
                reg_dt[i] = newgeno_dt
            master_dt['reg_dt'] = reg_dt     
            self._fn_lst = list(np.unique(self._fn_lst))
            for curfn in self._fn_lst:
                pd.DataFrame([json.dumps(master_dt)]).to_hdf(os.path.join(self.curdn, curfn), key='master_dt')
            
            if self.verbose: print('\nDone')

        def save_linkage_region(self, *, i, fn, keyfmt='ld/chrom{chrom}/i{i}/{varname}'): 
            # using 'store' instead of 'save' to indicate a connected relationship with 
            # the files used for this storage.
            geno_dt = self.reg_dt[i]
            chrom = geno_dt['chrom']
            curdn = self.curdn
            store_dt = dict() #geno_dt['store_dt']
            for varname, var in geno_dt.items():
                if varname in self._save_vars:
                    curfn  = fn.format(**locals())
                    key    = keyfmt.format(**locals())
                    var    = geno_dt[varname]
                    vartype = type(var)
                    if vartype is np.ndarray: vartype = var.dtype.type
                    curfullfn = os.path.join(curdn,curfn)
                    pd.DataFrame(var).to_hdf(curfullfn, key=key)
                    file_dt = dict(fn=curfn, key=key, 
                                   typestr=vartype.__module__+'.'+vartype.__name__)
                    store_dt[varname] = file_dt
                    self._fn_lst.append(curfn)
                    if self.verbose: print(f'saving: fn={curfullfn} key={key}'+' '*50,end='\r')
            geno_dt['store_dt'] = store_dt
                   
    # SumStat: ##############
    if True:

        def retrieve_sumstats_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_sumstats_region(i=i)

        def retrieve_sumstats_region(self, *, i):
            geno_dt = self.reg_dt[i] 
            sst_df  = geno_dt['sst_df']
            if 'beta_mrg' in geno_dt.keys():
                return None # Sumstat present so no need to compute anything.
            geno_dt['beta_mrg'] = self.compute_sumstats_region(i=i)
            if not 'beta_mrg' in sst_df.columns:
                geno_dt['sst_df']['beta_mrg'] = geno_dt['beta_mrg']
                
        retrieve_betamrg_region = retrieve_sumstats_region

        def retrieve_ldscores_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.retrieve_ldscores_region(i=i)

        def retrieve_ldscores_region(self, *, i):
            geno_dt = self.reg_dt[i]
            sst_df = geno_dt['sst_df']
            if not 'lds' in sst_df.columns:
                newsst_df = self.compute_ldscores_region(i=i)
                geno_dt['sst_df'] = newsst_df
            if self.clear_linkage:
                self.clear_linkage_region(i=i)

    # Clearing Functions: #####
    if True:

        def clear_all_xda(self):
            while len(self.xda_q) != 0:
                i_2_rm, key = self.xda_q.popleft()
                if i_2_rm == -1:
                    continue  # Continue to next iter if encountering a padding -1
                rmgeno_dt = self.reg_dt[i_2_rm]
                self.cur_total_size_in_gb -= getsizeof(rmgeno_dt[key].val) / 1024 ** 3
                rmgeno_dt.pop(key)
            [self.xda_q.append((-1,'')) for _ in range(5)]  # put 5x -1 in queue
            
        def clear_linkage_allregions(self):
            for i, geno_dt in self.reg_dt.items():
                self.clear_linkage_region(i=i)
            if self.verbose: print('\nDone')

        def clear_linkage_region(self, *, i):
            geno_dt = self.reg_dt[i]
            key_lst = list(geno_dt.keys())
            for key in key_lst:
                if key in self._clear_vars:
                    geno_dt.pop(key)
            if self.verbose: print(f'Cleared linkage region #{i} on chr{geno_dt["chrom"]}', end='\r'); sys.stdout.flush()
            

    ############################
    ## Get: ####################
    
    # Local Linkage: ###########
    if True:

        def get_auto_linkage_region(self, *, i):
            return self.get_specificied_linkage_region(i=i, shiftletter='D')

        def get_left_linkage_region(self, *, i):
            return self.get_specificied_linkage_region(i=i, shiftletter='L')

        def get_right_linkage_region(self, *, i):
            return self.get_specificied_linkage_region(i=i, shiftletter='R')

        def get_specificied_linkage_region(self, *, i, shiftletter):
            try:
                return self.reg_dt[i][shiftletter]
            except KeyError as e:
                if self._onthefly_retrieval:
                    if '_glocal' in shiftletter:
                        self.retrieve_linkage_region_glocalshiftwindow(i=i)
                    elif shiftletter in 'LDR':
                        self.retrieve_linkage_region(i=i)
                    elif shiftletter == 'Z':
                        self.retrieve_linkage_region_global(i=i)
                    else:
                        raise Exception(f'shiftletter={shiftletter}, on-the-fly retrieval not a valid option.')
                    try:
                        return self.reg_dt[i][shiftletter]
                    except Exception as e:
                        print('Failed, eventough on-the-fly retrieval was attempted')
                        raise e
                else:
                    raise Exception('on-the-fly retrieval blocked, set _onthefly_retrieval=True if desired')

        def get_auto_range_region(self, *, i):
            return self.reg_dt[i]['start_j'], self.reg_dt[i]['stop_j']

        def get_left_range_region(self, *, i):
            return self.reg_dt[i]['start_j_L'], self.reg_dt[i]['stop_j_L']

        def get_right_range_region(self, *, i):
            return self.reg_dt[i]['start_j_R'], self.reg_dt[i]['stop_j_R']

    # Sumstats: #################
    if True:
        
        def get_s(self):
            sst_df = self.get_sumstats_cur()
            try: 
                s = sst_df[['s']].values
                assert np.isnan(s).sum() == 0
                return s
            except:
                stansda = self.get_stansda()
                s = self.get_stansda().stats[:,[1]]
                self.s = s
                return s
            
        def get_sumstats_cur(self):
            sst_df_lst = []
            for i, geno_dt in self.reg_dt.items():
                sst_df = geno_dt['sst_df']
                sst_df_lst.append(sst_df)
            sst_df = pd.concat(sst_df_lst, axis=0)
            return sst_df

        def get_stansda(self, standardizer='unit'):
            if not standardizer=='unit':
                raise NotImplementedError('contact dev')
            
            if hasattr(self, 'stansda'):
                if type(self.stansda) is UnitTrained:
                    return self.stansda
                else:
                    raise NotImplementedError('contact dev')
                    
            standardizer_list = []
            for i, geno_dt in self.reg_dt.items():
                #(not 'stansda' in geno_dt.keys())
                if (not type(geno_dt['stansda']) is UnitTrained) & self._onthefly_retrieval:
                    self.retrieve_linkage_region(i=i)
                if type(geno_dt['stansda']) is UnitTrained:
                    standardizer_list.append(geno_dt['stansda'])
                else:
                    raise Exception('No standardizer detected. Compute this first. Contact dev if issue persists.')

            assert np.all([type(stan) is UnitTrained for stan in standardizer_list])            
            sid = np.concatenate([stan.sid for stan in standardizer_list])
            assert np.unique(sid).shape[0] == sid.shape[0]

            stats = np.concatenate([stan.stats for stan in standardizer_list])
            combined_unit_standardizer = UnitTrained(sid, stats)
            self.stansda = combined_unit_standardizer
            return combined_unit_standardizer


        def get_beta_marginal_full(self):
            beta_mrg_lst = []
            for i, geno_dt in self.reg_dt.items():
                beta_mrg_lst.append(geno_dt['beta_mrg'])
            beta_mrg_full = np.concatenate(beta_mrg_lst)
            return beta_mrg_full

        get_beta_marginal = get_beta_marginal_full
        
        def get_beta_marginal_region(self, *, i):
            return self.reg_dt[i]['beta_mrg']

    # Xda: ######################
    if True:
    
        def get_sda(self, *, i):
            geno_dt = self.reg_dt[i]
            if 'sda' in geno_dt.keys():
                return geno_dt['sda']
            else:
                if 'srd' in geno_dt.keys():
                    sda = geno_dt['srd'].read(dtype=self.dtype)
                    sda, stansda = sda.standardize(standardizer=geno_dt['stansda'], return_trained=True)
                    geno_dt['sda'] = sda
                    geno_dt['stansda'] = stansda
                    if self._save_s2sst:
                        geno_dt['sst_df']['s'] = stansda.stats[:,[1]]

                    if 'loaded_sda' in geno_dt.keys():
                        self.reloaded_xda_cnt += 1
                        if self.reloaded_xda_cnt in [5, 20, 100, 400]:
                            warnings.warn(
                                f'Reloaded sda for the {self.reloaded_xda_cnt}\'th time. This causes memory swapping,'
                                ' that might make the computation of linkage quite slow.'
                                'Probably because memory limits and/or linkage size.')
                    # Size determination and accounting:
                    geno_dt['loaded_sda']=True
                    self.cur_total_size_in_gb += getsizeof(sda.val) / 1024 ** 3
                    self.xda_q.append((i,'sda'))  # put respective i in queue.
                    while self.cur_total_size_in_gb > self.gb_size_limit:  # Keep removing till size is ok
                        i_2_rm, key = self.xda_q.popleft()
                        if i_2_rm == -1:
                            continue  # Continue to next iter if encountering a padding -1
                        rmgeno_dt = self.reg_dt[i_2_rm]
                        self.cur_total_size_in_gb -= getsizeof(rmgeno_dt[key].val) / 1024 ** 3
                        rmgeno_dt.pop(key)
                        if len(self.xda_q) <= 4:
                            raise Exception('The memory footprint of current settings is too high, '
                                            'reduce blocksize and/or correction windows or increase memory limits.')
                    return sda
                else:
                    raise Exception(f'No srd or sda found in region i={i}, this is not supposed to happen.')

        def get_pda(self):
            if not hasattr(self, 'pda'):
                pda = self.prd.read(dtype=self.dtype)
                pda, self.stanpda = pda.standardize(return_trained=True,
                                standardizer=self.pda_standardizer())
                self.pda = pda
            return self.pda
    

''
class LinkageData(BaseLinkageData):
    pass

def load_bimfam(base_fn, strip=True, bim=True, fam=True, fil_arr=None):
    if strip and (base_fn.split('.')[-1] in ('bim','fam','bed')): base_fn = '.'.join(base_fn.split('.')[:-1])
    bim_df = pd.read_csv(base_fn + '.bim', delim_whitespace=True, header=None, 
                         names=['chrom', 'snp', 'cm', 'pos', 'A1', 'A2']) if bim else None
    fam_df = pd.read_csv(base_fn + '.fam', delim_whitespace=True, header=None, 
                         names=['fid', 'iid', 'father', 'mother', 'gender', 'trait']) if fam else None
    
    if fil_arr is not None:
        bim_df = bim_df[bim_df.snp.isin(fil_arr)]
        bim_df = bim_df.reset_index(drop=True)
        
    return bim_df, fam_df

def load_linkagedata(fn):
    curfn = glob.glob(fn.format_map(defaultdict(lambda:'*')))[-1]
    master_dt = json.loads(pd.read_hdf(curfn, key='master_dt').loc[0,0])
    master_dt['curdn'] = os.path.dirname(curfn)
    linkdata = LinkageData(master_dt=master_dt)
    linkdata.load_linkage_allregions()
    return linkdata
