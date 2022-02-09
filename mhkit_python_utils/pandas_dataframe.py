
from datetime import datetime

import pandas as pd
import numpy as np


def timeseries_to_pandas(ts,ind,x):
    
    if x>1:
        ts=list(map(list,zip(*ts)))       
        df=pd.DataFrame(data=ts,index=ind)        
    else:        
        df=pd.DataFrame(data=ts,index=ind)
        
    return df.astype('float64')


def spectra_to_pandas(frequency,spectra,x,cols=None):
    if x>1:       
        ts=list(map(list,zip(*spectra)))       
        df=pd.DataFrame(data=ts,index=frequency)  
    else:
        df=pd.DataFrame(data=spectra,index=frequency)
        df.indexname='(Hz)'
        c_name=['PM']
    if cols is not None: 
        df.columns = cols
    return df.astype('float64')


def lis(li,app):
    li.append(app)
    return li


def datetime_index_to_ordinal(df):
    
    def to_ordinal_fraction(x):
        day = x.toordinal()
        dt = datetime.combine(x.date(), datetime.min.time(), x.tzinfo)
        fraction = (x.to_pydatetime() - dt).total_seconds() / 86400.
        return day + fraction + 366
    
    return np.array(list(map(to_ordinal_fraction, df.index)))
