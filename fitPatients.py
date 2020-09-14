# ------------------------ Imports ------------------------
import pandas as pd
import numpy as np
import scipy
import sys
import os
from tqdm import tqdm
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from itertools import product
import multiprocessing as mp
import re
from lmfit import minimize, Parameters

sys.path.append('./utils/')
import LotkaVolterraModel as lvm
import myUtils as utils
from fittingUtils import residual, PerturbParams, ComputeRSquared, LoadPatientData, PatientToOutcomeMap, \
    GetBestFit, GenerateFitSummaryDf, GenerateFitSummaryDf_AllPatients, PlotParameterDistribution_PatientCohort, \
    PlotFits

# Format plots
sns.set(style="white",
        font_scale=1.5,
        rc={'figure.figsize':(12,6)})
# ================================ Script setup ==========================================
dataDir = "./data/clinicalData/Bruchovsky_et_al/"
patientIdList = [int(re.findall(r'\d+',x)[0]) for x in os.listdir(dataDir)]
modelList = ["4params","3params_noCost","3params_noTurnover","2params_noCost_noTurnover"]
nFits = 25
perturbICs = True
runParallel = True
nProcesses = 24
outDir = './data/fits'

# Parameterise fitting algorithm
eps_data = 1.
solver_kws={'method':'DOP853', 'absErr':1.0e-8, 'relErr':1.0e-6, 'suppressOutputB':True}
optimiser_kws = {'method':'least_squares', 'xtol':1e-8, 'ftol':1e-8,
                 'nan_policy':'omit', 'verbose':0}
params = Parameters()
params.add('rS', value=0.027, min=1e-4, max=0.1, vary=False)
params.add('cost', min=0, max=1., vary=True)
params.add('rR', expr='(1-cost)*rS', vary=False)
params.add('turnover', min=0, max=1., vary=True)
params.add('dS', expr='turnover*rS', vary=False) # Constrain d<r
params.add('dR', expr='dS', vary=False) # Constrain dS = dR
params.add('dD', value=1.5, min=1, max=2, vary=False)
params.add('k', value=1, vary=False)
params.add('D', value=1, vary=False)
params.add('theta', value=1, vary=False)
params.add('DMax', value=1, vary=False)
params.add('n0', min=0.1, max=1., vary=True)
params.add('fR', min=1e-5, max=0.25, vary=True)
params.add('S0', expr='n0*(1-fR)', vary=False)
params.add('R0', expr='n0*fR', vary=False)

# ============================= Auxillary Functions ======================================
def FitModel(job):
    patientId, fitId, params, outDir = job['patientId'], job['fitId'], job['params'], job['outDir']
    dataDf = LoadPatientData(patientId,dataDir)
    currOutDir = os.path.join(outDir,"patient%d"%(patientId))
    if os.path.isfile(os.path.join(currOutDir,"fitObj_patient_%d_fit_%d.p"%(patientId, fitId))): return 0
    utils.mkdir(currOutDir)
    seed = int.from_bytes(os.urandom(4), byteorder='little')
    np.random.seed(seed)
    tmpModel = lvm.LotkaVolterraModel()
    if perturbICs: params = PerturbParams(params)
    try:
        fitObj = minimize(residual, params, args=(0, dataDf, eps_data, tmpModel, "PSA", solver_kws),**optimiser_kws)

        # Plot best fit
        myModel = lvm.LotkaVolterraModel()
        myModel.SetParams(**fitObj.params.valuesdict())
        myModel.Simulate(treatmentScheduleList=utils.ExtractTreatmentFromDf(dataDf),max_step=1,**solver_kws)
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        plt.plot(dataDf.Time, dataDf.PSA, linestyle='none', marker='x')
        myModel.Plot(plotPops=True, ymin=0, legendB=False, ax=ax)
        plt.savefig(os.path.join(currOutDir,"patient_%d_fit_%d.png" % (patientId, fitId)))
        plt.close()

        # Save fit
        fitObj.patientId = patientId
        fitObj.fitId = fitId
        fitObj.seed = seed
        fitObj.eps_data = eps_data
        fitObj.rSq = ComputeRSquared(fitObj,dataDf)
        pickle.dump(obj=fitObj, file=open(os.path.join(currOutDir,"fitObj_patient_%d_fit_%d.p"%(patientId, fitId)), "wb"))
    except:
        pass

# ================================= Main ================================================
pool = mp.Pool(processes=nProcesses,maxtasksperchild=1) if runParallel else None

jobList = []
for model in modelList:
    currParams = params.copy()
    for setting in model.split("_"):
        if setting == "noCost":
            currParams['cost'].set(value=0,vary=False)
        if setting == "noTurnover":
            currParams['turnover'].set(value=0,vary=False)
    currOutDir = os.path.join(outDir,model)
    for patientId,fitId in tqdm(product(patientIdList,range(nFits)),disable=runParallel):
        job = {'patientId':patientId, 'fitId':fitId, 'params':currParams, 'outDir':currOutDir}
        jobList.append(job)
        if not runParallel: FitModel(job)
if runParallel: list(tqdm(pool.imap(FitModel, jobList), total=len(jobList)))

# Analyse data
paramList = ["cost","turnover","n0","fR"]
for model in modelList:
    fitDir = os.path.join(outDir,model)
    dataToAnalyse = GenerateFitSummaryDf_AllPatients(patientIdList,fitDir,dataDir=dataDir)
    dataToAnalyse.to_csv(os.path.join(fitDir, "fitSummaryDf.csv"))

    # Plot all fits
    fig, axList = plt.subplots(6, 12, sharex=True, sharey=True, figsize=(30, 15))
    for i, patientId in enumerate(dataToAnalyse.PatientId.unique()):
        fit = GetBestFit(patientId, fitDir=fitDir)
        PlotFits(fits=[fit.fitId], dataDf=LoadPatientData(patientId,dataDir=dataDir), fitDir=fitDir,
                 solver_kws={**solver_kws,'max_step':1.},
                 ylim=2.5, titleStr="P%d;R2=%1.2f" % (patientId, fit.rSq), ax=axList.flatten()[i])
    plt.savefig(os.path.join(fitDir,"overview.pdf"))

    # Plot the parameters
    PlotParameterDistribution_PatientCohort(dataToAnalyse,paramList=paramList,fitDir=fitDir,
                                            outName=os.path.join(fitDir,"paramDist.pdf"))