'''Script to calculate bulk composition of our first data set'''

import CH3F_calc
import numpy as np
import csv
import pandas as pd


# Data array using mass 35
# In situ comes first in pairs, then adduct line-corrected
# feb16_data = np.array([
# [43.057, 0.025, 37.841, 37.607, 0.009, -99.9, -78.42, 0.73, -194.58, 0.83], # Lignin #1
# [38.669, 0.030, 35.510, 35.431, 0.007, 1.18, 1.06, 0.51, -111.52, 0.90], # Bamboo #1
# [18.237, 0.016, 23.415, 23.317, 0.006, -0.3, -0.98, 0.76, -92.44, 1.82], # Syringic acid #1
# [45.297, 0.017, 42.278, 42.169, 0.006, 43.19, 56.14, 0.76, -107.995, 1.712], #MeOH 1
# [17.646, 0.035, 17.532, 17.577, 0.009, 73.32, 64.12, 0.53, 6.51, 1.20]]) # CH3I-F 2
#
# # Data array using mass 15
# feb16_data = np.array([
# [43.057, 0.025, 37.841, 37.607, 0.009, -99.9, -78.42, 0.73, -190.31, 0.97], # Lignin #1
# [38.669, 0.030, 35.510, 35.431, 0.007, 1.18, 1.06, 0.51, -108.34, 0.69], # Bamboo #1
# [18.237, 0.016, 23.415, 23.317, 0.006, -0.3, -0.98, 0.76, -91.40, 0.63], # Syringic acid #1
# [45.297, 0.017, 42.278, 42.169, 0.006, 43.19, 56.14, 0.76, -111.22, 0.47], #MeOH 1
# [17.646, 0.035, 17.532, 17.577, 0.009, 73.32, 64.12, 0.53, 3.011, 0.962]]) # CH3I-F 2

# Data array using both, error weighted
feb16_data = np.array([
[43.063, 0.024, 37.841, 37.607, 0.009, -161.954, -131.224, 0.671, -192.78, 0.40], # Lignin #1
[38.669, 0.030, 35.510, 35.441, 0.007, -60.847, -58.197, 0.486, -109.50, 0.30], # Bamboo #1
[18.237, 0.016, 23.415, 23.317, 0.006, -56.263, -46.331, 0.425, -91.51, 0.35], # Syringic acid #1
[45.297, 0.017, 42.278, 42.169, 0.006, -32.629, -28.21, 0.376, -111.00, 0.20], #MeOH 1
[17.646, 0.035, 17.529, 17.564, 0.009, 36.774, 35.422, 0.467, 4.38, 0.56]]) # CH3I-F 2



bulk_comps_full = np.zeros((np.shape(feb16_data)[0],6))
bulk_comps_DFS = np.zeros((np.shape(feb16_data)[0],6))
bulk_comps_both = np.zeros((np.shape(feb16_data)[0],6))

for i in range(np.shape(feb16_data)[0]):
    bulk_comps_DFS_temp = CH3F_calc.MC_error_estimator(feb16_data[i,0], feb16_data[i,1],
    feb16_data[i,6], feb16_data[i,7], dDFS = feb16_data[i,8], dDFS_se = feb16_data[i,9])

    bulk_comps_both_temp = CH3F_calc.MC_error_estimator(feb16_data[i,0], feb16_data[i,1],
    feb16_data[i,6], feb16_data[i,7], d35_full = feb16_data[i,3], d35_full_se = feb16_data[i,4],
    dDFS = feb16_data[i,8], dDFS_se = feb16_data[i,9])

    bulk_comps_full_temp = CH3F_calc.MC_error_estimator(feb16_data[i,0], feb16_data[i,1],
    feb16_data[i,6], feb16_data[i,7], d35_full = feb16_data[i,3], d35_full_se = feb16_data[i,4])

    bulk_comps_full[i,:] = np.concatenate((np.mean(bulk_comps_full_temp,0),np.std(bulk_comps_full_temp,0)),0)
    bulk_comps_DFS[i,:] = np.concatenate((np.mean(bulk_comps_DFS_temp,0),np.std(bulk_comps_DFS_temp,0)),0)
    bulk_comps_both[i,:] = np.concatenate((np.mean(bulk_comps_both_temp,0),np.std(bulk_comps_both_temp,0)),0)

    print('number {0} complete'.format(str(i)))


fileName = 'Feb16_calculated_summary'

bcd = pd.DataFrame(bulk_comps_DFS)
bcb = pd.DataFrame(bulk_comps_both)
bcf = pd.DataFrame(bulk_comps_full)
bcf.to_csv('bcf.csv')
bcd.to_csv('bcd.csv')
bcb.to_csv('bcb.csv')
