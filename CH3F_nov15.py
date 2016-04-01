'''Script to calculate bulk composition of our first data set'''

import CH3F_calc
import numpy as np
import csv
import pandas as pd


# # Data array using mass 35
# # In situ comes first in pairs, then adduct line-corrected
# nov15_data = np.array([
# [16.817, 0.022, 16.568,16.490, 0.009, 66.53,53.26, 0.71, 14.29, 2.20], # CH3I #1, using next day's adduct lines
# [44.15, 0.026, 39.684, 39.552, 0.009, 28.75, 26.0306, 0.49,  -123.35, 1.43], # MeOH #1, No mass 16 DFS measurements for this
# [np.nan, np.nan, 16.386, 16.224, 0.007, 63.16, 57.63, 0.45, np.nan, np.nan], #
# [44.857, 0.018, 40.983, 40.853, 0.007, 31.65, 30.34, 0.51, -119.41, 0.71], # MeOH #2, only mass 16 DFS measurements
# [-7.736, 0.028, -0.688, -0.688, 0.008, 185.21, 175.16, 0.67, np.nan, np.nan],
# [20.64, 0.019, 27.3, 27.287, 0.007, 275.01, 256.46, 0.54, 222.37, 2.26 ],  # Vanillin #2
# [16.148, 0.02, 16.401, 16.388, 0.007, 58.7, 56.61, 0.40, 9.98, 3.26] ]) # MeOH #2,
#
# # Data array using mass 15
#
# nov15_data = np.array([
# [16.817, 0.022, 16.568,16.490, 0.009, 66.53,53.26, 0.71, 16.86, 0.55], # CH3I #1, using next day's adduct lines
# [44.15, 0.026, 39.684, 39.552, 0.009, 28.75, 26.0306, 0.49,  -123.35, 1.43], # MeOH #1, No mass 16 DFS measurements for this
# [np.nan, np.nan, 16.386, 16.224, 0.007, 63.16, 57.63, 0.45, np.nan, np.nan],
# [44.857, 0.018, 40.983, 40.853, 0.007, 31.65, 30.34, 0.51, -119.41, 0.71], # MeOH #2, only mass 16 DFS measurements
# [-7.736, 0.028, -0.688, -0.688, 0.008, 185.21, 175.16, 0.67, np.nan, np.nan],
# [20.64, 0.019, 27.3, 27.287, 0.007, 275.01, 256.46, 0.54, 232.95, 0.90 ],  # Vanillin #2
# [16.148, 0.02, 16.401, 16.388, 0.007, 58.7, 56.61, 0.40, 10.85, 1.09] ]) # MeOH #2,
#
#
# # Data array using both, error weighted
#
# nov15_data = np.array([
# [16.817, 0.022, 16.568,16.490, 0.009, 66.53,53.26, 0.71, 16.71, 0.29], # CH3I #1, using next day's adduct lines
# [44.15, 0.026, 39.684, 39.552, 0.009, 28.75, 26.0306, 0.49, -123.35, 1.43], # MeOH #1, No mass 16 DFS measurements for this
# [np.nan, np.nan, 16.386, 16.224, 0.007, 63.16, 57.63, 0.45, np.nan, np.nan],
# [44.857, 0.018, 40.983, 40.853, 0.007, 31.65, 30.34, 0.51, -119.41, 0.71], # MeOH #2, only mass 16 DFS measurements
# [-7.736, 0.028, -0.688, -0.688, 0.008, 185.21, 175.16, 0.67, np.nan, np.nan],
# [20.64, 0.019, 27.3, 27.287, 0.007, 275.01, 256.46, 0.54, 231.51, 0.70 ],   # Vanillin #2
# [16.148, 0.02, 16.401, 16.388, 0.007, 58.7, 56.61, 0.40, 10.76, 1.08] ]) # MeOH #2,

nov15_data = np.array([
[16.817, 0.022, 16.568,16.490, 0.009, 38.399,33.480, 0.69, 16.71, 0.29], # CH3I #1, using next day's adduct lines
[44.15, 0.026, 39.685, 39.553, 0.010, -26.015, -24.259, 0.447, -123.35, 1.43], # MeOH #1, No mass 16 DFS measurements for this
[np.nan, np.nan, 16.386, 16.224, 0.007, 63.16, 57.63, 0.45, np.nan, np.nan],
[44.857, 0.018, 40.892, 40.853, 0.007, -22.922, -23.035, 0.404, -119.41, 0.71], # MeOH #2, only mass 16 DFS measurements
[-7.736, 0.028, -0.688, -0.688, 0.008, 185.21, 175.16, 0.67, np.nan, np.nan],
[20.644, 0.019, 27.289, 27.276, 0.007, 225.706, 222.296, 0.487, 231.51, 0.70 ],   # Vanillin #2
[16.148, 0.02, 16.401, 16.388, 0.007, 36.822, 35.888, 0.395, 10.76, 1.08] ]) # CH3I number #3,


bulk_comps_full = np.zeros((np.shape(nov15_data)[0],6))
bulk_comps_DFS = np.zeros((np.shape(nov15_data)[0],6))
bulk_comps_both = np.zeros((np.shape(nov15_data)[0],6))

for i in range(np.shape(nov15_data)[0]):
    bulk_comps_DFS_temp = CH3F_calc.MC_error_estimator(nov15_data[i,0], nov15_data[i,1],
    nov15_data[i,6], nov15_data[i,7], dDFS = nov15_data[i,8], dDFS_se = nov15_data[i,9])

    bulk_comps_both_temp = CH3F_calc.MC_error_estimator(nov15_data[i,0], nov15_data[i,1],
    nov15_data[i,6], nov15_data[i,7], d35_full = nov15_data[i,3], d35_full_se = nov15_data[i,4],
    dDFS = nov15_data[i,8], dDFS_se = nov15_data[i,9])

    bulk_comps_full_temp = CH3F_calc.MC_error_estimator(nov15_data[i,0], nov15_data[i,1],
    nov15_data[i,6], nov15_data[i,7], d35_full = nov15_data[i,3], d35_full_se = nov15_data[i,4])

    bulk_comps_full[i,:] = np.concatenate((np.mean(bulk_comps_full_temp,0),np.std(bulk_comps_full_temp,0)),0)
    bulk_comps_DFS[i,:] = np.concatenate((np.mean(bulk_comps_DFS_temp,0),np.std(bulk_comps_DFS_temp,0)),0)
    bulk_comps_both[i,:] = np.concatenate((np.mean(bulk_comps_both_temp,0),np.std(bulk_comps_both_temp,0)),0)

    print('number {0} complete'.format(str(i)))


fileName = 'nov15_calculated_summary'

bcd = pd.DataFrame(bulk_comps_DFS)
bcb = pd.DataFrame(bulk_comps_both)
bcf = pd.DataFrame(bulk_comps_full)
bcf.to_csv('bcf.csv')
bcd.to_csv('bcd.csv')
bcb.to_csv('bcb.csv')
