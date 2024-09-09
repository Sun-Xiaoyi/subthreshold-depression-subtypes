# -*- coding: utf-8 -*-

import os
import pcntoolkit as pcn

wdir1 = 'D:\\Data_Chen\\With_DIDA_all_HC\\data_forNorm\\'
covfile = os.path.join(wdir1, 'age_sex_HC_Final.txt')
respfile = os.path.join(wdir1, 'zFCS_HC_forNorm.txt')

testcov = os.path.join(wdir1, 'age_sex_StD_Final.txt')
testresp = os.path.join(wdir1, 'zFCS_StD_forNorm.txt')

# To assess model generalizability, we initially estimated normative models 
# in the healthy group using 10-fold cross-validation.
pcn.normative.estimate(covfile, respfile, cvfolds = 10, alg = 'gpr', outputsuffix = '_10fold')

# Estimated the normative models of FCS as a function of age and sex using GPR and 
# calculated individual FCS deviations in normative models for StD participants.
pcn.normative.estimate(covfile, respfile, testresp=testresp, testcov=testcov, alg="gpr", saveoutput='True')

