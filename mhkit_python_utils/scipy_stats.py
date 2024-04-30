

def convert_to_array(obj, x, method):
    '''
    Function used to convert scipy.stats.rv_continous object into
    dtypes that can be used in matlab

    Parameters
    ----------
    obj: scipy.stats.rv_continuous
        Scipy statistical object (frozen).
    x: np.array
        Array to analyze
    method: str
        Method to extract from the scipy statistical object. Possible
        choices include: pdf, cdf, ppf, sf
    
    Returns
    -------
    stat: array-like
        Array-like result from the scipy object based on the chosen method
    '''

    if method=='pdf':
        return obj.pdf(x)
    if method=='cdf':
        return obj.cdf(x)
    if method=='ppf':
        return obj.ppf(x)
    if method=='sf':
        return obj.sf(x)
    else:
        print('ERROR: no valid method selected!!!')
        raise