# -*- coding: utf-8 -*-
"""
Created on Thu May  5 10:05:30 2022

@author: goyal080
"""

import numpy as np
from scipy.stats import uniform
from scipy.stats import norm
import gurobipy as gp
from gurobipy import GRB
import os
                                
def mycallback(model, where):
    if where == GRB.Callback.MIPSOL:
        xval = model.cbGetSolution(model._x.values())
        S = [j for j in N if round(xval[j])==1]        
        NminusS = [j for j in N if round(xval[j])==0]        
        for i in range(m): 
            add_weak_cuts(model._x, model._w, S, NminusS, i, 1, [])
    if where == GRB.Callback.MIP:  
        crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)  
        if crt_cutcnt != model._cutcnt: 
            model._cutcnt = crt_cutcnt            
            
def mycallback2(model, where):
    if where == GRB.Callback.MIPSOL:
        xval = model.cbGetSolution(model._x.values())        
        S = [j for j in N if round(xval[j])==1]        
        NminusS = [j for j in N if round(xval[j])==0]        
        for i in range(m): 
            add_lifted_cuts1(model._x, model._w, S, NminusS, i, 1, True)
            add_lifted_cuts2(model._x, model._w, S, NminusS, i, 1, True)
    if where == GRB.Callback.MIP:  
        crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)  
        if crt_cutcnt != model._cutcnt: 
            model._cutcnt = crt_cutcnt 

def mycallback3(model, where):
    if where == GRB.Callback.MIPSOL:
        xval = model.cbGetSolution(model._x.values())
        S = [j for j in N if round(xval[j])==1]
        NminusS = [j for j in N if round(xval[j])==0]
        for i in range(m): 
            add_lifted_cuts1(model._x, model._w, S, NminusS, i, 1, False)
            add_lifted_cuts2(model._x, model._w, S, NminusS, i, 1, False)                    
    if where == GRB.Callback.MIP:  
        crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)  
        if crt_cutcnt != model._cutcnt:             
            model._cutcnt = crt_cutcnt                   

#%
def add_weak_cuts(x, w, S, NminusS, i, lazy, cnt):
    a_Set = lambda Set: sum(a[j,i] for j in Set)
    a_S = a_Set(S)
    a_N = a_Set(N)
    
    if lazy == 1:
        model.cbLazy(w[i] <= f(a_S) - sum((f(a_N)-f(a_N-a[j,i]))*(1-x[j]) for j in S) + sum((f(a_S+a[j,i])-f(a_S))*x[j] for j in NminusS))
        model.cbLazy(w[i] <= f(a_S) - sum((f(a_S)-f(a_S-a[j,i]))*(1-x[j]) for j in S) + sum((f(a[j,i])-f(0))*x[j] for j in NminusS))
    else:
        model.addConstr(w[i] <= f(a_S) - sum((f(a_N)-f(a_N-a[j,i]))*(1-x[j]) for j in S) + sum((f(a_S+a[j,i])-f(a_S))*x[j] for j in NminusS), name='ineq1_'+str(i)+'_'+str(cnt))
        model.addConstr(w[i] <= f(a_S) - sum((f(a_S)-f(a_S-a[j,i]))*(1-x[j]) for j in S) + sum((f(a[j,i])-f(0))*x[j] for j in NminusS), name='ineq2_'+str(i)+'_'+str(cnt))

#%
def find_interval(A_lb, A_ub, delta):
    test_ub = A_ub-delta >= 0
    test_lb = A_lb-delta <= 0
    result = np.where(test_ub&test_lb==True)[0]
    if result.size > 0:
        return result[0]
    else:
        return None

#%
def add_lifted_cuts1(x, w, S, NminusS, i, lazy, approx):     
    a_Set = lambda Set: sum(a[j,i] for j in Set)

    def zeta(delta):
        if delta <= A1_lb[-1]:
            k = sz_NminusS-1
        else:
            k = find_interval(A1_lb, A1_ub, delta)
        func_val = f(a_S+A[k]+delta) - sum((f(a_S+a_sort[j])-f(a_S)) for j in range(k+1)) - f(a_S)
        return func_val  
    
    def gamma(delta):
        if delta <= A2_lb[-1]:
            return zeta(delta)
        else:
            k = find_interval(A2_lb, A2_ub, delta)
            if k == None or k==0:
                return zeta(delta)
            else:
                func_val = zeta(A2_ub[k]) - (f(a_S+a_sort[k])-f(a_S))*(A2_ub[k]-delta)/a_sort[k]
                return func_val
            
    a_S = a_Set(S)
    a_N = a_Set(N)        
    
    if S == [] or NminusS == []:
        if lazy == 1:
            model.cbLazy(w[i] <= f(a_S) - sum((f(a_N)-f(a_N-a[j,i]))*(1-x[j]) for j in S) + sum((f(a_S+a[j,i])-f(a_S))*x[j] for j in NminusS))
    else:
        a_sort = -np.sort(-a[NminusS, i])
        sz_NminusS = len(a_sort)
        A = np.cumsum(a_sort)        
        # %--------------------------------------------------------------------%
        A1_lb = -A
        A1_ub = np.insert(-A[0:sz_NminusS-1], 0, 0)               
        mu = []
        for k in range(sz_NminusS):
            mu.append(g(gprime_inv(a_sort[k]/(f(a_S+a_sort[k])-f(a_S))))-a_S)
        A2_lb = mu-A
        A2_ub = mu-np.insert(A[0:sz_NminusS-1], 0, 0)       
        # %--------------------------------------------------------------------%
        if lazy == 1:
            if approx == True:
                model.cbLazy(w[i] <= f(a_S) + sum(gamma(-a[j,i])*(1-x[j]) for j in S) + sum((f(a_S+a[j,i])-f(a_S))*x[j] for j in NminusS))
            else:
                model.cbLazy(w[i] <= f(a_S) + sum(zeta(-a[j,i])*(1-x[j]) for j in S) + sum((f(a_S+a[j,i])-f(a_S))*x[j] for j in NminusS))                

#%
def add_lifted_cuts2(x, w, S, NminusS, i, lazy, approx):   
    a_Set = lambda Set: sum(a[j,i] for j in Set)
        
    def ksi(delta):
        if delta >= A3_ub[-1]:
            k = sz_S-1
        else:
            k = find_interval(A3_lb, A3_ub, delta)
        func_val = f(a_S-A[k]+delta) + sum((f(a_S)-f(a_S-a_sort[j])) for j in range(k+1)) - f(a_S)
        return func_val   
    
    def omega(delta):
        if delta >= A4_ub[-1]:
            return ksi(delta)
        else:
            k = find_interval(A4_lb, A4_ub, delta)
            if k == None or k==0:
                return ksi(delta)
            else:
                func_val = ksi(A4_ub[k]) - (f(a_S)-f(a_S-a_sort[k]))*(A4_ub[k]-delta)/a_sort[k]
                return func_val 
    
    a_S = a_Set(S)
    a_N = a_Set(N)
            
    if S == [] or NminusS == []:
        if lazy == 1:
            model.cbLazy(w[i] <= f(a_S) - sum((f(a_S)-f(a_S-a[j,i]))*(1-x[j]) for j in S) + sum((f(a[j,i])-f(0))*x[j] for j in NminusS))
    else:
        a_sort = -np.sort(-a[S, i])
        sz_S = len(a_sort)
        A = np.cumsum(a_sort)        
        # %--------------------------------------------------------------------%
        A3_lb = np.insert(A[0:sz_S-1], 0, 0)  
        A3_ub = A              
        nu = []
        for k in range(sz_S):
            nu.append(a_S-g(gprime_inv(a_sort[k]/(f(a_S)-f(a_S-a_sort[k])))))
        A4_lb = np.insert(A[0:sz_S-1], 0, 0) - nu
        A4_ub = A -nu       
        # %--------------------------------------------------------------------% 
        if lazy == 1:
            if approx == True:
                model.cbLazy(w[i] <= f(a_S) - sum((f(a_S)-f(a_S-a[j,i]))*(1-x[j]) for j in S) + sum(omega(a[j,i])*x[j] for j in NminusS))  
            else:
                model.cbLazy(w[i] <= f(a_S) - sum((f(a_S)-f(a_S-a[j,i]))*(1-x[j]) for j in S) + sum(ksi(a[j,i])*x[j] for j in NminusS))                   
     
#%
def optimize_model(n, m1, lmbd, instance, Time_limit, exact_lifting):
    
    global N, a, m, model, f, g, gprime_inv
    
    m = m1
    N = list(range(n))
    
    np.random.seed(instance*100)
    a_constr = uniform.rvs(0, 0.2, size = n)
    alpha = uniform.rvs(0.05, 0.10, size = n)
    beta = uniform.rvs(0, 1, size = n)
    r = np.zeros([n,m])
    a = np.zeros([n,m])
    prob_pi = np.ones(m)/m
    
    for i in range(m):
        ln_f = norm.rvs(0.05, np.sqrt(0.0025))  # normal distribution
        eps = norm.rvs(0, np.sqrt(0.0025), size = n)  # normal distribution
        r[:,i] = np.exp(alpha + beta*ln_f + eps)
        a[:,i] = np.multiply(r[:,i], a_constr)
       
    f = lambda y: -np.exp(-y/lmbd)
    g = lambda z: -lmbd*np.log(-z)
    gprime_inv = lambda y: -lmbd/y
        
    #%   
    filename = "LogFiles/Weak"+"_"+str(n)+"_"+str(m)+"_"+str(lmbd)+"_inst"+str(instance)+".log"
    if os.path.exists(filename):    os.remove(filename)
    env = gp.Env()
    env.start()
    model = gp.Model('cap_budget1', env)
    x = model.addVars(n, vtype=GRB.BINARY, name="x")
    w = model.addVars(m, vtype=GRB.CONTINUOUS, lb=-float('inf'), ub=0, name="w")
    model._x = x
    model._w = w 
    model._cutcnt = 0
    model.addConstr(sum(a_constr[j]*x[j] for j in N) <= 1, name="budget_constr")
    model.setObjective(1 + sum(prob_pi[i]*w[i] for i in range(m)), GRB.MAXIMIZE)
    model.Params.Cuts = 0
    model.Params.LazyConstraints = 1
    model.Params.TimeLimit = Time_limit   
    model.Params.LogToConsole = 0
    model.Params.LogFile = filename
    model.optimize(mycallback)   
    model.Params.LogFile = ""
    env.close()
    return_vec = {}
    return_vec[1] = [model._cutcnt, model.NodeCount, model.ObjVal, model.ObjBound, model.MIPGap*100, model.Runtime]
            
    #%
    filename = "LogFiles/Lifted_Approx_"+"_"+str(n)+"_"+str(m)+"_"+str(lmbd)+"_inst"+str(instance)+".log"
    if os.path.exists(filename):    os.remove(filename)
    env = gp.Env()
    env.start()
    model = gp.Model('cap_budget2', env)
    x = model.addVars(n, vtype=GRB.BINARY, name="x")
    w = model.addVars(m, vtype=GRB.CONTINUOUS, lb=-float('inf'), ub=0, name="w")
    model._x = x
    model._w = w 
    model._cutcnt = 0
    model.addConstr(sum(a_constr[j]*x[j] for j in N) <= 1, name="budget_constr")
    model.setObjective(1 + sum(prob_pi[i]*w[i] for i in range(m)), GRB.MAXIMIZE)
    model.Params.Cuts = 0
    model.Params.LazyConstraints = 1
    model.Params.TimeLimit = Time_limit   
    model.Params.LogToConsole = 0
    model.Params.LogFile = filename
    model.optimize(mycallback2)   
    model.Params.LogFile = ""
    env.close()
    return_vec[2] = [model._cutcnt, model.NodeCount, model.ObjVal, model.ObjBound, model.MIPGap*100, model.Runtime]  
    
    #%
    if exact_lifting == True:
        filename = "LogFiles/Lifted_Exact_"+"_"+str(n)+"_"+str(m)+"_"+str(lmbd)+"_inst"+str(instance)+".log"
        if os.path.exists(filename):    os.remove(filename)
        env = gp.Env()
        env.start()
        model = gp.Model('cap_budget3', env)
        x = model.addVars(n, vtype=GRB.BINARY, name="x")
        w = model.addVars(m, vtype=GRB.CONTINUOUS, lb=-float('inf'), ub=0, name="w")
        model._x = x
        model._w = w 
        model._cutcnt = 0
        model.addConstr(sum(a_constr[j]*x[j] for j in N) <= 1, name="budget_constr")
        model.setObjective(1 + sum(prob_pi[i]*w[i] for i in range(m)), GRB.MAXIMIZE)
        model.Params.Cuts = 0
        model.Params.LazyConstraints = 1
        model.Params.TimeLimit = Time_limit   
        model.Params.LogToConsole = 0
        model.Params.LogFile = filename
        model.optimize(mycallback3)   
        model.Params.LogFile = ""
        env.close()
        return_vec[3] = [model._cutcnt, model.NodeCount, model.ObjVal, model.ObjBound, model.MIPGap*100, model.Runtime]  
    
    return return_vec
