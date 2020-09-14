# ================================================================================================================
# Script to obtain TTP for CT compared to AT on a cohort of tumours with different initial densities and sensitive
# fractions.
# ================================================================================================================
import numpy as np
import pandas as pd
import sys
from itertools import product
from tqdm import tqdm
import multiprocessing as mp
# Import my own libraries
sys.path.append('./utils/')
from odeAnalysisUtils import rdModel_nsKill,GenerateParameterDic,Simulate_ContinousTx,ProfileTreatmentStrategies

# Parameterise script
intervalLength = 1.
atThreshold = 0.5
turnoverList = np.linspace(0., 0.5, 251)
costList = np.linspace(0.,0.5, 251)
nProcesses = 20
# ------------------------ Functions ------------------------
def ProfileParamSet(paramSet):
    turnover,cost = paramSet
    paramDic = {"rS": .027, "rR": .027, "cRS": 1., "cSR": 1., "dD": 1.5,
                "k": 1., "D": 0, "theta": 1, 'DMax': 1.}
    # Set up the parameters
    _, currParamDic = GenerateParameterDic(initialSize=0, rFrac=0,
                                           cost=cost, turnover=turnover,
                                           paramDic=paramDic)
    print(cost,turnover)
    # Perform the comparison
    txComparisonDf = ProfileTreatmentStrategies(modelFun=rdModel_nsKill,paramDic=currParamDic,
                                                enableProgressBar=False,
                                                atThresholdList=[atThreshold],intervalLength=intervalLength,
                                                initialSizeList=[0.25],
                                                rFracList=[0.001],
                                                tumourSizeWhenProgressed=1.2,cureThreshold=1e-7)

    # Save results
    txComparisonDf['Turnover'] = turnover
    txComparisonDf['Cost'] = cost
    return txComparisonDf

# ------------------------ Main ------------------------
pool = mp.Pool(processes=nProcesses,maxtasksperchild=1)
jobList = list(product(turnoverList, costList))
tmpDicList = list(tqdm(pool.imap(ProfileParamSet,jobList),total=len(jobList)))
resultsDf = pd.concat(tmpDicList,sort=True)
resultsDf.to_csv("ttpAnalysis_costTurnoverSweep.csv")