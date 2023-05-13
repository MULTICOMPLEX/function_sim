from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import math

from tqdm import tqdm

# Reference: Cover & Thomas "Elements of Information Theory" pp. 333-336
# Reference: Justin Dauwels "Numerical Computation of the Capacity of Continuous Memoryless Channels"
#            http://www.dauwels.com/files/memoryless.pdf
# Erik Leitinger Capacity and Capacity-Achieving Input Distribution of the Energy Detector https://www.spsc.tugraz.at/sites/default/files/master_thesis_ED.pdf
# Chan, Hranilovic, Kschischang: Capacity-Achieving Probability Measure for Conditionally Gaussian Channels With Bounded Inputs https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1435651
# Joel Smith (1972) The Information Capacity of Amplitude and Variance-Constrained Scalar Gaussian Channels / https://ac.els-cdn.com/S0019995871903469/1-s2.0-S0019995871903469-main.pdf?_tid=aefc5895-6825-4e58-929e-41d8e1a82bb1&acdnat=1552018847_d07922b271c821875314e04643825997

# discretization = .02

R_min = 0
R_max = 40
eps = 1e-7
check_every = 1
# get_capacity is very slow so we do lots of
# optimization and not much checking where we are.

# this is a lot faster than scipy
# ... and apparently inlining is faster than this.
def normal_cdf(x, mean, std):
  return 0.5 * (1 + math.erf((x - mean) / (std * 1.4142135623730951)))

def get_std(mean):
  return 0.0037 * mean ** 2 + 0.0247 * mean + 0.0437

# channel has x (the conditioning variable) in rows and y in columns. 
def get_channel(discretization, increments):
  channel = np.zeros((increments, increments))
  for ix, i in enumerate(tqdm(np.arange(R_min, R_max, discretization))):
    mean = i + discretization/2
    std = get_std(mean)
    for jx, j in enumerate(np.arange(R_min, R_max, discretization)):
      xl = j
      xu = j + discretization
      # apologies for the inlining for speedup
      # density = norm.cdf(xu, loc=mean, scale=std) - norm.cdf(xl, loc=mean, scale=std)
      density = 0.5 * ((1 + math.erf((xu - mean) / (std * 1.4142135623730951))) - (1 + math.erf((xl - mean) / (std * 1.4142135623730951))))
      channel[ix][jx] = density
    # channel[ix] = channel[ix] / np.sum(channel[ix]) # since the distributions are truncated
    # print ((1 - np.sum(channel[ix])) / increments)
    if (1 - np.sum(channel[ix])) > 1e-6:
      channel[ix] += (1 - np.sum(channel[ix])) / increments
  return channel

# TODO vectorize more
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

def get_capacity(q, r, channel):
  # caps = np.zeros_like(channel)
  capacity = 0
  for x in range(q.shape[1]):
    for y in range(q.shape[0]):
      if r[x] != 0 and channel[x][y] != 0:
        # caps[x][y] = r[x] * channel[x][y] * np.nan_to_num((np.log2(q.T)[x][y] - np.log2(r)[x]))
        capacity += r[x] * channel[x][y] * (np.log2(q[y][x]) - np.log2(r[x]))
  # print (caps)
  # print (r * channel * np.nan_to_num((np.log2(q.T) - np.log2(r))))
  return capacity
  # return np.sum(r * channel * np.nan_to_num((np.log2(q.T) - np.log2(r))))

def plot(channel, vals=[2.25,3,4.25,6,8.5,13,25]):
  plt.figure()
  x_axis = np.arange(channel.shape[0])*discretization + R_min
  for r_mean in vals:
    bucket = int((r_mean - R_min) / discretization)
    plt.plot(x_axis, channel[bucket])
  plt.show()

def BA_discrete(discretization):
  increments = int((R_max - R_min)/discretization)
  channel = get_channel(discretization, increments)
  # plot(channel)

  capacity = 0
  # q has y (the conditioning variable) as rows and x as columns.
  q = np.ones((increments, increments))
  q = q/q.shape[0]
  r = np.ones((increments))
  r = r/r.size
  i = 0
  while True:
    i += 1
    # print ("q",q)
    # print ("r",r)
    q = update_q(q, r, channel)
    r = update_r(q, r, channel)
    if i % check_every == 0:
      new_capacity = get_capacity(q, r, channel)
      print (new_capacity)
      yield np.arange(R_min, R_max, discretization), r
      # if new_capacity - capacity < eps * check_every:
      #   return np.arange(R_min, R_max, discretization), r
      capacity = new_capacity

# fig, ax = plt.subplots()
# iterator = BA_discrete(0.2)
# x1, y1 = next(iterator)
# line, = ax.plot(x1, y1)


def get_next(i):
  x, y = next(iterator)
  line.set_ydata(y)
  global check_every
  if check_every < 10:
    check_every += 1
  else:
    check_every = round(check_every*1.1)
  print (check_every)
  return line

def main():
  # print (BA_discrete(0.002))
  # res = [BA_discrete(10*2**(-i)) for i in range(10)]
  # x, y = BA_discrete(0.8)
  # plt.plot(x, y)

  anim = FuncAnimation(fig, get_next, frames=np.arange(0, 105), interval=100)
  fig.suptitle('Blahut-Arimoto Optimization')
  plt.xlabel('Resistance (kOhm)')
  plt.ylabel('Probability')
  anim.save('dist_0.2.gif', dpi=80, writer='imagemagick')

  # res = [1.694730274542409, 2.3095430902786878, 2.8401561849782375, 3.3037163404128953, 3.690441508981898, 3.962557410364222, 4.0776002325659775, 4.1449106180324415, 4.178185139995351, 4.188598885819961]
  # print (res)
  # plt.plot([40/(10*2**(-i)) for i in range(10)], res)
  plt.xscale("log")

  plt.savefig("dist_1.png")

def main2():
  res = [1.694730274542409, 2.3095430902786878, 2.8401561849782375, 3.3037163404128953, 3.690441508981898, 3.962557410364222, 4.0776002325659775, 4.1449106180324415, 4.178185139995351, 4.188598885819961]
  # plt.plot(x, y)
  plt.plot([40/(10*2**(-i)) for i in range(10)], res, label="empirical capacity")
  plt.plot([40/(10*2**(-i)) for i in range(10)], [3 for _ in range(10)], label="current 8-level model")
  plt.title('Blahut-Arimoto Optimization')
  plt.xlabel('Number of discretization segments')
  plt.ylabel('Capacity (bits)')
  plt.xscale("log")
  plt.legend()
  plt.savefig("discr.png")

if __name__ == "__main__":
  main2()