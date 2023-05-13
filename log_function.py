#Monte Carlo Log function 
import numpy as np

def mc_log(t, n):
    x = np.random.uniform(1, t, n)
    total = np.sum(1/x)
    return (t - 1) * total / n

# Example usage:
t = .5
n = 10000000
approx_ln = mc_log(t, n)
print(f"ln({t}) = {approx_ln:.6f}")
print("exact =", np.log(t))