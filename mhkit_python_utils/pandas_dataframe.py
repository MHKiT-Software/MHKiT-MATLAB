import pandas as pd


def timeseries_to_pandas(ts,ind,x):
    
    if x>1:
        ts=list(map(list,zip(*ts)))       
        df=pd.DataFrame(data=ts,index=ind)        
    else:        
        df=pd.DataFrame(data=ts,index=ind)
        
    return df.astype('float64')

def spectra_to_pandas(frequency,spectra,x):
    if x>1:       
        ts=list(map(list,zip(*spectra)))       
        df=pd.DataFrame(data=ts,index=frequency)  
    else:
        df=pd.DataFrame(data=spectra,index=frequency)
        df.indexname='(Hz)'
        c_name=['PM']
    return df.astype('float64')

def lis(li,app):
    li.append(app)
    return li
    
    

