# -*- coding: utf-8 -*-

from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio
coord_file = "D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\res\\vol_coor.txt"
output_dir = "D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\res"

filenames = volume(coord_file, output_dir)

brain_map = "D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\gene_abagen\\res\\vox_z_clus1.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], resample=True)
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_z_clus1.mat',{'surrogate_maps':surrogate_maps})
