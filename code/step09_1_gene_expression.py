# -*- coding: utf-8 -*-

import abagen
import pandas as pd
import numpy as np

expression,report = abagen.get_expression_data('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\code\\L_shen268_group.nii',
                                               atlas_info=pd.read_excel('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\code\\L_shen268_group_info.xlsx'),
                                               return_donors=True,
                                               ibf_threshold=0.5,
                                               data_dir='D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\data\\microarray\\',
                                               return_report=True)

expression_filter,stability = abagen.keep_stable_genes(list(expression.values()),
                                                             threshold=0.1,percentile=False,return_stability=True)
expression_filter = pd.concat(expression_filter).groupby('label').mean()
np.savetxt('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\res\\L_stability.csv', stability, delimiter=',')
expression_filter.to_csv('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\res\\L_expression_filter.csv')
