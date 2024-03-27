# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 14:00:22 2023

@author: pwilson3

Updated and sorted list of all the folders with PPC data. Can be called to
other scripts.
"""

# The location of the folders containing the raw measurement data.
directory = "C:\\Users\\louis\\OneDrive - University of Ottawa\\uOttawa\\2023-24\\Hiver2024\\Sunab pt 2\\Data\\"

# Place the names of the folders that contain the raw .sciv files in these 
# lists.
# You can create as many lists as you want. The list objects will be callable 
# from other files once you import this script.
# The data_collector file is 
# 10J measurements
"""
folder_C5200 = ['20230316_C5200_X19Y16', '20230322_C5200_X21Y16', 
                '20230322_C5200_X22Y14', '20230323_C5200_X21Y7', 
                '20230329_C5200_X22Y14', '20230329_C5200_X21Y13',
                '20230330_C5200_X21Y13', '20230908_C5200_X21Y14_Gavin',
                '20230920_C5200_X19Y13', '20231120_C5200_X19Y13',
                '20231221_C5200_X26Y13', '20231221_C5200_X27Y10'] 
"""


# 2J measurements
folder_C5195 = ['20240221_C5195_X25Y5','20240221_C5195_X19Y15']

folder_C5245_25 = ['20240109_C5245_X3Y1','20240118_C5245_X6Y1', '20240326_C5245_X7Y0-25degIV'] #'20240130_C5245_X7Y0-25' ,'20240306_C5245_X7Y0'

folder_C5245_27 = ['20240314_C5245_X7Y0-27degIV']#, '20240326_C5245_X7Y0-27degFixedCurrent'

fixed_current_23 = ['20240326_C5245_X7Y0-23degFixedCurrent']

fixed_current_27 = ['20240326_C5245_X7Y0-27degFixedCurrent']

folder_C5245_23 = ['20240306_C5245_X7Y0-23degIV_wrong_parameters_in_files']#,'20240326_C5245_X7Y0-23degFixedCurrent'

folder_C5246 = ['20240124_C5246_X12Y1', '20240111_C5246_X8Y1']

folder_C5247 = ['20240116_C5247_X5Y4', '20240116_C5247_X7Y5','20240123_C5247_X5Y5', '20240123_C5247_X7Y6']
