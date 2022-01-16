import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from osgeo import gdal, osr
from gdalconst import *


def biweeklyDatasets(path, isoflag, chemstryFlag):

    """
        Parameters:
            - airtemp: (array of) daily average air temperatures [Celsius].
            - rh: (array of) daily average relative humidity [%].
            - airpress: (array of) daily average air pressure data [Pa].
            - u: (array of) daily average wind speed at 2 m [m/s].
            - Rnet: (array of) net radiation [J/m2/d]
    """

    biweekly_df = pd.read_excel(path + 'peatland_data.xlsx', sheet_name='Data', skiprows=[1,2]).iloc[18:,1:]
    timeIndexPlusOne = np.unique(biweekly_df['Datum'].to_numpy())[:-1]
    timeIndexPlusOne = timeIndexPlusOne[np.where(np.logical_not(np.isnan(timeIndexPlusOne)))]
    timeIndex = timeIndexPlusOne[1:]

    if isoflag=='2H':
        isoKeyWord = '2H'
    elif isoflag=='18O':
        isoKeyWord = '18O'

    # ********** Biweekly data extraction *************
    delta_sw = np.full((len(timeIndex), 9), np.nan)
    no3_sw = np.full((len(timeIndex), 9), np.nan)
    nh4_sw = np.full((len(timeIndex), 9), np.nan)
    srp_sw = np.full((len(timeIndex), 9), np.nan)
    doc_sw = np.full((len(timeIndex), 9), np.nan)
    cl_sw = np.full((len(timeIndex), 9), np.nan)
    fe_sw = np.full((len(timeIndex), 9), np.nan)
    ca_sw = np.full((len(timeIndex), 9), np.nan)
    delta_gw = np.full((len(timeIndex), 6), np.nan)
    no3_gw = np.full((len(timeIndex), 6), np.nan)
    nh4_gw = np.full((len(timeIndex), 6), np.nan)
    srp_gw = np.full((len(timeIndex), 6), np.nan)
    doc_gw = np.full((len(timeIndex), 6), np.nan)
    cl_gw = np.full((len(timeIndex), 6), np.nan)
    fe_gw = np.full((len(timeIndex), 6), np.nan)
    ca_gw = np.full((len(timeIndex), 6), np.nan)



    if chemstryFlag==False:
        for i in range(len(timeIndex)):
            df_sw = biweekly_df.iloc[i*18:(i+1)*18][:9]
            df_gw = biweekly_df.iloc[i*18:(i+1)*18][10:16]
            delta_sw[i, :] = df_sw[isoKeyWord]
            delta_gw[i, :] = df_gw[isoKeyWord]

    elif chemstryFlag:
        for i in range(len(timeIndex)):
            df_sw = biweekly_df.iloc[i*18:(i+1)*18][:9]
            df_gw = biweekly_df.iloc[i*18:(i+1)*18][10:16]
            delta_sw[i, :] = df_sw[isoKeyWord]
            no3_sw[i, :] = df_sw['NO3--N']
            nh4_sw[i, :] = df_sw['NH4+-N']
            srp_sw[i, :] = df_sw['SRP']
            doc_sw[i, :] = df_sw['DOC']
            cl_sw [i, :] = df_sw['Cl-']
            fe_sw [i, :] = df_sw['Fe']
            ca_sw [i, :] = df_sw['Ca']
            delta_gw[i, :] = df_gw[isoKeyWord]
            no3_gw[i, :] = df_gw['NO3--N']
            nh4_gw[i, :] = df_gw['NH4+-N']
            srp_gw[i, :] = df_gw['SRP']
            doc_gw[i, :] = df_gw['DOC']
            cl_gw [i, :] = df_gw['Cl-']
            fe_gw [i, :] = df_gw['Fe']
            ca_gw [i, :] = df_gw['Ca']


    return timeIndexPlusOne, delta_sw, no3_sw, nh4_sw, srp_sw, doc_sw, cl_sw, fe_sw, ca_sw, delta_gw, no3_gw, nh4_gw, srp_gw, doc_gw, cl_gw, fe_gw, ca_gw

def dailyDatasets(path, isoflag, timeIndex):
    _df = pd.read_excel(path+'climate_daily.xlsx', index_col='TIMESTAMP')
    P = np.full(len(timeIndex)-1, np.nan)
    deltaP = np.full(len(timeIndex) - 1, np.nan)
    E0 = np.full(len(timeIndex)-1, np.nan)
    for i in range(len(timeIndex)-1):
        df = _df[timeIndex[i]:timeIndex[i+1]][1:]
        dailyP = df['Rain_corr_mm_Tot']
        dailyDeltaP = df[isoflag+'Rain']
        dailyE0 = df['E0']
        dailyE0[dailyE0<0]=0
        P[i] = np.sum(dailyP)
        deltaP[i] = np.sum(dailyP*dailyDeltaP)/P[i]
        E0[i] = np.sum(dailyE0)

    return P, deltaP, E0

def IDW_deltaGW(path, delta_gw, delta_sw):

    from inverseDistanceWeightedInterpolation import tree

    df = pd.read_excel(path + 'peatland_data.xlsx', sheet_name='Info')
    corr_gw = df.iloc[1:7, 8:10].to_numpy()
    corr_sw = df.iloc[1:10, 3:5].to_numpy()
    deltaRin = np.full(np.shape(delta_sw), np.nan)

    spacingx = np.linspace(446975, 447975, 6)
    spacingy = np.linspace(5806025, 5808025, 12)
    X2 = np.meshgrid(spacingx, spacingy)
    grid_shape = X2[0].shape
    X2 = np.reshape(X2, (2, -1)).T
    corrIdx_sw = np.full(np.shape(corr_sw), -9999)
    corrIdx_sw[:,0] = np.around(((corr_sw[:,0]- 446975)/200).astype(float))
    corrIdx_sw[:,1] = np.around(((corr_sw[:,1]-5806025)/200).astype(float))

    delta_gw[[8,10,12],[0,4,5]] = np.nan
    for j in range(np.shape(delta_gw)[1]):
        idx_notNan = np.where(np.logical_not(np.isnan(delta_gw[:,j])))
        f = interpolate.interp1d(idx_notNan[0], delta_gw[:,j][idx_notNan], kind='linear', fill_value="extrapolate")
        idx_isNan = np.where(np.isnan(delta_gw[:,j]))
        delta_gw[:,j][idx_isNan] = f(idx_isNan[0])


    for i in range(np.shape(delta_gw)[0]):
        nanIdx = np.isnan(delta_gw[i,:])
        if nanIdx.all():
            continue
        X = corr_gw[np.where(nanIdx==False)]
        Z = delta_gw[i,:][np.where(nanIdx==False)]
        idw_tree = tree(X, Z, leafsize=10)


        z2 = idw_tree(X2)
        z3 = z2.reshape(grid_shape)
        deltaRin[i,:] = z3[corrIdx_sw[:,1],corrIdx_sw[:,0]]

    return deltaRin

def GIS(path):

    #***************** Sort the cumulative length along the stream in 1 m raster **************
    ds = gdal.Open(path + 'GIS/mainStreamRaster.tif', GA_ReadOnly)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    arr = arr.astype(np.double)

    arr[arr==15] = np.nan
    arr[arr<15] = 1
    arr[0,39:42] = np.nan
    tmp = np.argwhere(np.nansum(arr, axis=1)==2)

    for i in tmp[:,0]:

        idx = np.argwhere(arr[i,:]==1)[:,0]
        arr[i, idx[0]] = np.nan
        arr[i, idx[1]] = 1.4

    tmp = np.nansum(arr, axis=1)

    for j in range(np.shape(arr)[0]):
        arr[j,np.where(arr[j,:]>0)] = np.sum(tmp[:j])

    X = np.nansum(arr, axis=1)

    # ********** Interpolate width and depth along the stream ********
    corr_streamDepthWidth_project = np.loadtxt(path + 'GIS/streamDepthWidth.txt', skiprows=1, delimiter=',')
    xlist = corr_streamDepthWidth_project[:,-4]
    ylist = corr_streamDepthWidth_project[:,-3]
    depthList = corr_streamDepthWidth_project[:,6]/100
    widthList = corr_streamDepthWidth_project[:,7]
    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    prosrs = ds.GetSpatialRef()
    prosrs.ImportFromWkt(ds.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    ycorr_streamDepthWidth = np.array([])
    for i in range(len(xlist)):
        coords = ct.TransformPoint(xlist[i], ylist[i])
        xOffset = int((coords[0] - xOrigin) / pixelWidth)
        yOffset = int((coords[1] - yOrigin) / pixelHeight)
        ycorr_streamDepthWidth = np.append(ycorr_streamDepthWidth, yOffset)


    ycorr_streamDepthWidth[ycorr_streamDepthWidth<0] = np.nan
    ycorr_streamDepthWidth[ycorr_streamDepthWidth>np.shape(arr)[0]] = np.nan
    notNan_idx = np.where(np.logical_not(np.isnan(ycorr_streamDepthWidth)))

    X_samples = X[ycorr_streamDepthWidth[notNan_idx].astype(int)]


    f_depth = interpolate.interp1d(X_samples, depthList[notNan_idx], kind='linear', fill_value="extrapolate")
    Y_depth = f_depth(X)
    f_width = interpolate.interp1d(X_samples, widthList[notNan_idx], kind='linear', fill_value="extrapolate")
    Y_width = f_width(X)

    QGauge = np.array([[447510.645, 5808131.804], [447709.109, 5806275.576]])
    tmp1 = np.array([39.87755208, 76.67445833])  # water level depth at PN and PS on 04.05.2021
    tmp2 = np.array([36.31370833, 72.089625])  # water level depth at PN and PS on 26.04.2021
    tmp3 = np.array([54.957981, 49.827571])  # DEM at PN and PS from drone image on 26.04.2021
    ycorr_QGauge = np.array([])
    for j in range(2):
        xOffset = int((QGauge[j,0] - xOrigin) / pixelWidth)
        yOffset = int((QGauge[j,1] - yOrigin) / pixelHeight)
        ycorr_QGauge = np.append(ycorr_QGauge, yOffset)
    f_waterSurface = interpolate.interp1d(X[ycorr_QGauge.astype(int)], tmp3-tmp2+tmp1, kind='linear', fill_value="extrapolate")
    Y_waterSurface = f_waterSurface(X)

    Z_channelBottom = Y_waterSurface - Y_depth
    Y_channelBottom = Y_width - Y_depth


    # ************ Reach Length *************
    df = pd.read_excel(path + 'peatland_data.xlsx', sheet_name='Info')
    corr_sw = df.iloc[1:10, 3:5].to_numpy()
    ycorr_sw = ((corr_sw[:,1] - yOrigin) / pixelHeight).astype(int)

    reach_length = X[ycorr_sw][1:] - X[ycorr_sw][:-1]
    reachNumber = 8
    reach_channelBottom = np.full(reachNumber, np.nan)
    reach_depth = np.full(reachNumber, np.nan)
    reach_area = np.full(reachNumber, np.nan)
    reach_storage = np.full(reachNumber, np.nan)
    for k in range(reachNumber):
        x1 = int(X[ycorr_sw][k])
        x2 = int(X[ycorr_sw][k+1])
        reach_channelBottom[k] = np.sum(Y_channelBottom[x1:x2]*tmp[x1:x2])/np.sum(tmp[x1:x2])
        reach_depth[k] = np.sum(Y_depth[x1:x2]*tmp[x1:x2])/np.sum(tmp[x1:x2])
        reach_area[k] = np.sum((Y_channelBottom[x1:x2] + Y_depth[x1:x2])*tmp[x1:x2])
        reach_storage[k] = np.sum((Y_channelBottom[x1:x2]*2 + Y_depth[x1:x2]) * Y_depth[x1:x2] * 0.5 * tmp[x1:x2])

    print(reach_channelBottom, reach_depth, reach_area, reach_storage)

    """
    df_output1 = pd.DataFrame([])
    df_output2 = pd.DataFrame([])
    df_output3 = pd.DataFrame([])
    df_output1['X'] = X
    df_output1['depth'] = Y_depth
    df_output1['width'] = Y_width

    df_output2['Xsamples'] = X_samples
    df_output2['Zsamples'] = depthList[notNan_idx]
    df_output2['Ysamples'] = widthList[notNan_idx]

    df_output3['X_sw'] = X[ycorr_sw]
    df_output3['Y_sw'] = f_width(X[ycorr_sw])
    df_output3['Z_sw'] = f_depth(X[ycorr_sw])

    df_output1.to_csv('output1.csv')
    df_output2.to_csv('output2.csv')
    df_output3.to_csv('output3.csv')
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(6, 6))
    ax1.plot(X, Y_width)
    ax1.plot(X, Y_channelBottom)
    ax1.scatter(X_samples, widthList[notNan_idx])
    ax2.plot(X, Y_waterSurface)
    ax2.plot(X, Z_channelBottom)
    ax2.scatter(X_samples, f_waterSurface(X_samples) - depthList[notNan_idx])
    plt.show()

    """

    return Y_channelBottom, Z_channelBottom

def dataProcessing(path, isoFlag, chemstryFlag, plotFlag):

    timeIndexPlusOne, delta_sw, no3_sw, nh4_sw, srp_sw, doc_sw, cl_sw, fe_sw, ca_sw, delta_gw, no3_gw, nh4_gw, srp_gw,\
    doc_gw, cl_gw, fe_gw, ca_gw = biweeklyDatasets(path, isoFlag, chemstryFlag)

    P, deltaP, E0 = dailyDatasets(path, isoFlag, timeIndexPlusOne)

    deltaRin = IDW_deltaGW(path, delta_gw, delta_sw)

    if plotFlag:
        fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, sharey=True, figsize=(10,6))
        im1 =ax1.imshow(deltaP.reshape(1,-1), vmin=-70, vmax=-40)
        im2 = ax2.imshow(delta_sw.T, vmin=-70, vmax=-40)
        im3 = ax3.imshow(deltaRin.T, vmin=-70, vmax=-40)
        plt.show()

    upstream_delta_sw = np.full(np.shape(delta_sw), np.nan)
    upstream_delta_sw[:, 1:] = delta_sw[:, :-1]
    checkArr = np.full(np.shape(delta_sw), 1)
    for i in range(np.shape(delta_sw)[0]-1):
        idx1 = np.where(delta_sw[i,:]>deltaP[i])
        checkArr[i, idx1] = 0

        idx2 = np.where(delta_sw[i,:]>upstream_delta_sw[i,:])
        checkArr[i, idx2] = 0

        idx3 = np.where(delta_sw[i,:]>deltaRin[i,:])
        checkArr[i, idx3] = 0

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, sharey=True, figsize=(10, 6))
    im1 = ax1.imshow(deltaP.reshape(1, -1), vmin=-70, vmax=-40)
    im2 = ax2.imshow(delta_sw.T, vmin=-70, vmax=-40)
    im3 = ax3.imshow(deltaRin.T, vmin=-70, vmax=-40)

    im4 = ax4.imshow(checkArr.T)
    plt.show()







path = r'D:\OneDrive\Phd\WorkPackages\3\data/'
isoFlag = '2H'
plotFlag = False
chemstryFlag = False

#dataProcessing(path, isoFlag, chemstryFlag, plotFlag)
GIS(path)
