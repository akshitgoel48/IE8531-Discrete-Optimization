# -*- coding: utf-8 -*-
"""
Created on Sat May  7 15:40:27 2022

@author: goyal080
"""

import pandas as pd
import os
from capital_budgeting import optimize_model

if not os.path.exists('LogFiles'):
    os.mkdir('LogFiles')
      
column_list = [ '', '', '', '',
               '', '', 'Submodular', 'Ineq', '', '', '', 
               '', '', 'Approx', 'Lifiting', '', '', '',
               '', '', 'Exact', 'Lifitng', '', '']    

column_list2 = [ 'n', 'm', 'lambda', 'inst',
               '# cuts', '# nodes', 'ObjVal', 'ObjBnd', '%Gap', 'Runtime', '', 
               '# cuts', '# nodes', 'ObjVal', 'ObjBnd', '%Gap', 'Runtime', '',
               '# cuts', '# nodes', 'ObjVal', 'ObjBnd', '%Gap', 'Runtime']

tb_array = []
tb_array.append(column_list2)

TimeLimit = 600
exact_lifting = False

for n in [15]: # [15, 25, 50]:
    for m in [100]: #[1, 25, 100]:
        for lmbd in [1]: #[1, 2]:
            for instance in range(5):
                return_vec = optimize_model(n, m, lmbd, instance+1, TimeLimit, exact_lifting)
                if instance == 0:
                    tb_arr = [n, m, lmbd, instance+1]
                else:
                    tb_arr = ['', '', '', instance+1]                                        
                tb_arr.extend(return_vec[1])
                tb_arr.extend([''])
                tb_arr.extend(return_vec[2])
                if exact_lifting == True:
                    tb_arr.extend([''])
                    tb_arr.extend(return_vec[3])
                tb_array.append(tb_arr)
            tb_array.append(['']*len(column_list))
            
            table = pd.DataFrame(tb_array, columns = column_list) 
            table.to_excel('final_results.xlsx', index = False)
