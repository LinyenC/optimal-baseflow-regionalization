# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 18:19:48 2025

@author: 11519
"""


import pandas as pd
import numpy as np
import joblib
from sklearn.ensemble import RandomForestRegressor,RandomForestClassifier


df = pd.read_excel(r'.\stationInfo.xlsx',dtype = {'stationID':str})
allX = df[['Area', 'RE', 'vf_clay_s', 'Ksat', 'Wsat', 'TMAXmax', 'P_frac', 'P5', 'P50', 'P95', 'PEmax', 'RS_mean', 'RS_std']].values

notNan_index = ~np.isnan(allX).any(1)

RF2 = joblib.load(r'.\RF2.pkl')
RF3 = joblib.load(r'.\RF3.pkl')
RF4 = joblib.load(r'.\RF4.pkl')

ymodel_RF2 = RF2.predict(allX[notNan_index])
ymodel_RF3 = RF3.predict(allX[notNan_index])
ymodel_RF3 = np.round(ymodel_RF3,0)
ymodel_RF4 = RF4.predict(allX[notNan_index])
ymodel_RF4 = np.round(ymodel_RF4,0)
ymodel_RF4[ymodel_RF2==0] = ymodel_RF3[ymodel_RF2==0]

ymodel = np.full([df.shape[0]], np.nan)
ymodel[notNan_index] = ymodel_RF4

outdf = pd.DataFrame(ymodel, columns=['pred_N'], index = df.stationID)
outdf.to_excel(r'.\predicted_N.xlsx')