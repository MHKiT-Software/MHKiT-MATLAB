
from datetime import datetime

import pandas as pd
import numpy as np
import scipy.io
import pdb
def timeseries_to_pandas(ts,ind,x):
    
    if x>1:
        ts=list(map(list,zip(*ts)))       
        df=pd.DataFrame(data=ts,index=ind)        
    else:        
        df=pd.DataFrame(data=ts,index=ind)
        
    return df.astype('float64')

def convert_number_array_and_index_to_dataframe(columns_as_lists, index, transpose_threshold=1):
    """
    Convert a time series into a Pandas DataFrame.

    Parameters:
    time_series (list of lists): The time series data represented as a list of lists.
    index (list or array-like): The index or datetime index for the DataFrame.
    transpose_threshold (int): A threshold that controls whether to transpose the time series.

    Returns:
    pd.DataFrame: A Pandas DataFrame containing the time series data.

    If the 'transpose_threshold' is greater than 1, the time series is transposed before creating
    the DataFrame using the provided index. If the 'transpose_threshold' is not greater than 1,
    the DataFrame is created from the original time series without transposing, using the
    provided index. The data in the DataFrame is treated as float64 for consistency.

    Example:
    time_series = [[1, 2, 3], [4, 5, 6]]
    index = ['2023-01-01', '2023-01-02']
    df = convert_timeseries_to_dataframe(time_series, index, transpose_threshold=2)
    """

    if transpose_threshold > 1:
        time_series = list(map(list, zip(*time_series)))
        df = pd.DataFrame(data=time_series, index=index)
    else:
        df = pd.DataFrame(data=time_series, index=index)

    return df.astype('float64')


def list_to_series(input_list, index=None):
    return pd.Series(input_list, index)

def spectra_to_pandas(frequency,spectra,x,cols=None):
    if x>1:       
        ts=list(map(list,zip(*spectra)))  
        print("ts = ", ts)     
        df=pd.DataFrame(data=ts,index=frequency)  
        print("df = ", df)     

    else:
        df=pd.DataFrame(data=spectra,index=frequency)
        df.indexname='(Hz)'
        c_name=['PM']
    if cols is not None: 
        df.columns = cols
        print("df.astype('float64') = ", df.astype('float64'))
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

