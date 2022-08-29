#/Users/johannawahn/Desktop/iGEM/Modeling
# -*- coding: utf-8 -*-
"""
Created on %(29/6/2020)s

@author: %(johananWahn)s
"""
#%% PRIME NUMBER TEST
#%% FUNCTIONS
def is_prime(nbr):
    if (nbr == 1):
        return("is not prime")
    
    if (nbr == 2):
        return("is prime")
    
    for n in range(2,nbr):
        if nbr%n ==  0  :
            return("is not prime")
            break
        
        if n == (nbr-1):
            return("is prime")
        


#%%USER INPUT
txt = input("Type number you want to test: ")
print(txt, is_prime(int(txt)))


#%%FIBONACCI SEQUENCE Xn= Xn-1 + Xn-2
#%%FUNCTION
def fibonacci(nbr):
    seq = list()
    if (nbr == 0 or nbr == 1 ):
        return(list(range(nbr+1)))
    
    else:
        n1 = 0
        seq.append(n1)
        n2 = 1
        seq.append(n2)
        
        for n in range(0,nbr-1):
            Fn = n1 + n2
            seq.append(Fn)
            
            n1=n2
            n2=Fn
        
    return(seq)
            
            
sequence = fibonacci(7)            
    

    
    
    
    
    
    
    
    