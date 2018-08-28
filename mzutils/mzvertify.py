

def date_vertify(a):

    """" Take an array of dates, and remove/replace things that are not dates with a very old date 
    This really help to avoid date comparisons from toppling over """

    import pandas as pd
    import datetime

    b=[]

    for i in a:
            j=datetime.datetime(1900,1,1)
            try:
                    j = pd.DatetimeIndex([i])[0]
            except ValueError as e:
                    j = pd.DatetimeIndex([j])[0]
                    pass
            except KeyError as e:
                    j = pd.DatetimeIndex([j])[0]
                    pass
            except TypeError as e:
                    j = pd.DatetimeIndex([j])[0]
                    pass
            b.append(j)

            #print(j)

    return b



def date_compare(f):

    """ Compares dates in a dataframe; comares date in column 2 to col1, col 3 to col 2 and so on, and make sure each
    date following is higher than the previous one. If it is not, it warns and prints, but continues with unaltered data frame """


    for c in range(0,(f.shape[1]-1)):
        for r in range(1,f.shape[0]):
            if (type(f.iloc[r,c])==type(f.iloc[r,(c+1)])):
                if f.iloc[r,c+1] < f.iloc[r,c]:
                    print(["\nWarning these dates may be incorrect\n",f.iloc[r,],"\n\n"])




#def fillna_downbet(df):
#    import numpy as np
#    df = df.copy()
#    for col in df:
#        non_nans = df[col][~df[col].apply(np.isnan)]
#        start, end = non_nans.index[0], non_nans.index[-1]
#        df[col].loc[start:end] = df[col].loc[start:end].fillna(method='ffill')
#    return df


