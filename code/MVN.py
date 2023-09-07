"""
Compute the feasibility domain size using a multivariate normal
"""

import scipy 
from scipy.stats import multivariate_normal
import numpy as np

def feasibility(A): 
  V = np.linalg.inv(A)
  S=len(V[0,:])
  bounds=np.zeros(S)
  cov=100*np.dot(V,V.T)  #a factor 100 to spread out the distribution (somehow helps). 
  mean=np.zeros(S)
  return 2**S*multivariate_normal.cdf(bounds,mean,cov, allow_singular=True)
 
