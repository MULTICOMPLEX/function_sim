from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
import math

from tqdm import tqdm

np.random.seed(7)

R_min = 0
R_max = 40

KL_descent_count = 20

sqrt_2_pi = math.sqrt(2*math.pi)

def q(x):
  return 0.0037 * x * x + 0.0247 * x + 0.0437

def diff_q(x):
  return 0.0074 * x + 0.0247

pdf_cache = {}
def channel_pdf(y, x):  # p(y|x) on the channel model pdf
  if (x,y) in pdf_cache:
    return pdf_cache[(x,y)]
  qx = q(x)
  exp_frac = ((y-x)/2*qx)
  result = np.exp(-(exp_frac*exp_frac))/(sqrt_2_pi*qx)
  pdf_cache[(x,y)] = result
  if len(pdf_cache) > 50000:
    pdf_cache.clear()
  return result

def channel_pdf_diff_wrt_x(y, x):
  qx = q(x)
  dqx = diff_q(x)
  exp_frac = ((y-x)/2*qx)
  exp_term = math.exp(-(exp_frac*exp_frac))
  ft_num = dqx*(y-x)**2/(qx**3) - (y-x)/(qx*qx)
  st_num = qx*dqx*exp_term
  return ((exp_term*ft_num) - st_num/(qx*qx)/(sqrt_2_pi*qx))

# TODO this is some sketchy math that is worth re-deriving
def D_KL_diff_wrt_xi(y_evals, p_k_y, xi, wi):
  s = 0
  y_eval, H = y_evals
  for yl, pky in zip(y_eval, p_k_y):
    diff_x = channel_pdf_diff_wrt_x(yl, xi)
    if diff_x < 1e-50:
      continue
    py_x = channel_pdf(yl, xi)
    mult = 1 + math.log(py_x / pky) - wi * (py_x / pky)
    s += H * diff_x * mult
  return s

def initialize_particles(num_particles):
  L = [np.random.random()*(R_max - R_min) + R_min for _ in range(num_particles)]
  w = [1/len(L) for _ in range(num_particles)]
  return L, w

# channel has x (the conditioning variable) in rows and y in columns. 
def get_channel_matrix(L, pdf):
  channel_matrix = np.zeros((len(L), len(L)))
  for i, Li in enumerate(L):
    for j, Lj in enumerate(L):
      channel_matrix[i, j] = pdf(Lj, Li)
    channel_matrix[i] = channel_matrix[i] / np.sum(channel_matrix[i])
  return channel_matrix


def update_q(q, r, channel):
  q_new = np.zeros_like(q)
  for y in range(q_new.shape[0]):
    q_new[y] = (r * channel[:, y])
    q_new[y] /= np.sum(q_new[y])
  return q_new

def update_r(q, r, channel):
  r_new = np.prod(np.power(q.T, channel), axis=1)
  r_new /= np.sum(r_new)
  return r_new

# TODO vectorize?
def get_capacity(q, r, channel):
  capacity = 0
  for x in range(q.shape[1]):
    for y in range(q.shape[0]):
      if r[x] != 0 and channel[x][y] != 0:
        capacity += r[x] * channel[x][y] * (np.log2(q[y][x]) - np.log2(r[x]))
  return capacity

def BA_discrete(channel, L, r_init, eps=1e-8):
  capacity = 0
  check_every = 1
  r = r_init
  # q has y (the conditioning variable) as rows and x as columns.
  q = np.ones((len(r), len(r)))
  q = q/q.shape[0]
  i = 0
  while True:
    i += 1
    q = update_q(q, r, channel)
    r = update_r(q, r, channel)
    if i % check_every == 0:
      new_capacity = get_capacity(q, r, channel)
      # print (new_capacity)
      if new_capacity - capacity < eps * check_every:
        break
      capacity = new_capacity
  return r, capacity

def check_convergence(curr_cap, old_cap, eps=1e-7):
  return abs(curr_cap-old_cap) < eps

def KL_div_p_y_xi_pky(xi, p_k_y, pdf, y_evals):
  y_eval, H = y_evals
  s = 0
  for yl, pky in zip(y_eval, p_k_y):
    py_x = pdf(yl, xi)
    if py_x != 0:
      s += H * py_x * math.log(py_x / pky)
  return s

# TODO backtracking search
def update_once(xi, wi, p_k_y, y_evals, pdf, alpha=0.00175, beta=0.6):
  old_value = KL_div_p_y_xi_pky(xi, p_k_y, pdf, y_evals)
  gradient = D_KL_diff_wrt_xi(y_evals, p_k_y, xi, wi)
  for i in range(100):  # just to cap max number of iterations
    xi_new = xi + alpha * gradient
    value = KL_div_p_y_xi_pky(xi_new, p_k_y, pdf, y_evals)
    # print ("i", i, "xi", xi, "gradient", gradient, "value", value, "old_value", old_value, "alpha", alpha)
    if value > old_value and xi_new >= R_min and xi_new <= R_max:
      return xi_new, wi
    alpha = alpha * beta
  return xi, wi

def get_p_k_y(L, w, pdf, y_evals):
  y_eval, _ = y_evals
  p_k_y = np.zeros_like(y_eval)
  for i, y in enumerate(y_eval):
    p_k_y[i] = sum(wi * pdf(y, xi) for xi, wi in zip(L, w))
  p_k_y = p_k_y / np.sum(p_k_y)
  return p_k_y

def run_continuous_BA(num_particles, y_evals):
  # step 1
  L, w = initialize_particles(num_particles)
  curr_cap = 0
  while True:
    particle_channel_matrix = get_channel_matrix(L, channel_pdf)
    old_cap = curr_cap
    # step 2
    w, curr_cap = BA_discrete(particle_channel_matrix, L, w)
    if check_convergence(curr_cap, old_cap):
      break
    print (curr_cap)
    # step 3 (now below)
    # p_k_y = get_p_k_y(L, w, y_evals)
    # step 4
    for _ in tqdm(range(KL_descent_count)):
      for i, xi in enumerate(L):
        p_k_y = get_p_k_y(L, w, channel_pdf, y_evals)
        L[i], w[i] = update_once(xi, w[i], p_k_y, y_evals, channel_pdf)
  return curr_cap
  

def main():
  H = 0.04
  y_evals = (np.arange(R_min+H, R_max, H), H)
  print (run_continuous_BA(50, y_evals))

if __name__ == "__main__":
  main()
