import math
import decimal
from decimal import Decimal as D
from alive_progress import alive_bar
from contextlib import nullcontext
from scipy import stats
import multiprocessing
import tqdm
import functools
import monteflib

decimal.getcontext().prec = 100

FREQUENCIES = {
  "A": 0.0812,
  "B": 0.0149,
  "C": 0.0271,
  "D": 0.0432,
  "E": 0.1202,
  "F": 0.0230,
  "G": 0.0203,
  "H": 0.0592,
  "I": 0.0731,
  "J": 0.0010,
  "K": 0.0069,
  "L": 0.0398,
  "M": 0.0261,
  "N": 0.0695,
  "O": 0.0768,
  "P": 0.0182,
  "Q": 0.0011,
  "R": 0.0602,
  "S": 0.0628,
  "T": 0.0910,
  "U": 0.0288,
  "V": 0.0111,
  "W": 0.0209,
  "X": 0.0017,
  "Y": 0.0211
}

PREFERENCES = ["E", "T"]

ORDERED_LETTERS = PREFERENCES + [l for l in FREQUENCIES.keys() if l not in PREFERENCES]
ORDERED_FREQUENCIES = [FREQUENCIES[l] for l in ORDERED_LETTERS]
FREQSIZE = len(ORDERED_FREQUENCIES)
PREFSIZE = len(PREFERENCES)
PFINAL = D(1 - sum(ORDERED_FREQUENCIES))

trial_cache = {}
def load_trial_cache():
  try:
    with open("monte_carlo_cache", "r") as f:
      while l := f.readline():
        csv = [int(e) for e in l.replace("\n", "").split(",") if e.isdigit()]
        trial_cache[csv[0]] = [csv[1] / csv[2], csv[1], csv[2]]
  except FileNotFoundError:
    pass
load_trial_cache()

def save_trial_cache():
  sorted_keys = sorted(trial_cache.keys())
  with open("monte_carlo_cache", "w+") as f:
    for k in sorted_keys:
      f.write(",".join((str(k), str(trial_cache[k][1]), str(trial_cache[k][2]), str(trial_cache[k][0]))) + "\n")

prob_cache = {}
def load_cache():
  try:
    with open("cache", "r") as f:
      while l := f.readline():
        splitline = l.split("-")
        prob_cache[int(splitline[0])] = D(splitline[1])
  except FileNotFoundError:
    pass
load_cache()

def save_cache():
  sorted_keys = sorted(prob_cache.keys())
  with open("cache", "w+") as f:
    for k in sorted_keys:
      f.write(str(k) + "-" + str(prob_cache[k]) + "\n")

def get_computed_prob(length):
  if length in prob_cache:
    return prob_cache[length]
  prob = compute_prob(length)
  prob_cache[length] = prob
  save_cache()
  return prob

def compute_prob(length):
  return factorial(length) * recursive_sum(length)

sum_cache = {}

def recursive_sum(length, index=0, counts=[]):
  if index < PREFSIZE:
    prefsremaining = PREFSIZE - index - 1
    mintodistributeleftovers = math.floor((length + prefsremaining * (FREQSIZE + (1 - prefsremaining) / 2)) / (FREQSIZE + 1)) + 1
    lowbound = max(1 + prefsremaining, mintodistributeleftovers)
    # At minimum, remaining n prefs need 1 + 2 + 3 + ... + n (assuming strict inequality)
    highbound = length - (prefsremaining) * (prefsremaining + 1) // 2
    if counts:
      highbound = min(highbound, counts[-1] - 1)
  else:
    lettersafter = FREQSIZE - index
    lowestpref = counts[PREFSIZE - 1]
    lowbound = max(0, length - lettersafter * (lowestpref - 1))
    highbound = min(length, lowestpref - 1)

    cache_key = (lowestpref, length, index)
    if cache_key in sum_cache:
      return sum_cache[cache_key]

  tot = 0

  if index == 0:
    bar = alive_bar(highbound - lowbound + 1)
  else:
    bar = nullcontext()
  with bar as bar:
    for c in range(lowbound, highbound + 1):
      nextiter = recursive_sum(length - c, index + 1, counts + [c]) if index + 1 < FREQSIZE else 1
      toadd = nextiter * exponent(index, c) / factorial(c)
      if index + 1 == FREQSIZE:
        toadd *= exponent(-1, length - c) / factorial(length - c)
      tot += toadd
      if index == 0:
        bar()
  if index >= PREFSIZE:
    sum_cache[cache_key] = tot
  return tot

exponent_cache = [[D(1)] for i in range(len(ORDERED_FREQUENCIES) + 1)]
def exponent(index, power):
  if power >= len(exponent_cache[index]):
    if index == -1:
      prob = PFINAL
    else:
      prob = D(ORDERED_FREQUENCIES[index])
    for i in range(len(exponent_cache[index]), power + 1):
      exponent_cache[index].append(exponent_cache[index][i - 1] * prob)
  return exponent_cache[index][power]

factorial_cache = [1]
def factorial(x):
  if x >= len(factorial_cache):
    for i in range(len(factorial_cache), x + 1):
      factorial_cache.append(i * factorial_cache[-1])
  return factorial_cache[x]

def fortran_single_process(length, trialsize):
  return monteflib.trial(length, trialsize)

def brute_probability(length, trials, processnum=8, batchsize=10000):
  total_trials = trials // batchsize * batchsize

  total_success = 0
  batches = total_trials // batchsize
  process_func = functools.partial(fortran_single_process, length)
  with multiprocessing.Pool(processes=processnum) as pool:
    for r in tqdm.tqdm(pool.imap_unordered(process_func, [batchsize] * batches), total = batches):
      total_success += r
      pass
  return total_success / total_trials, total_success, total_trials

def smart_search(target_p):
  target = D(target_p)
  
  first_pos = -1
  L = 1
  while first_pos < 0:
    L *= 2
    print("Searching %s" % (L))
    if get_computed_prob(L) >= target:
      first_pos = L
      break

  return binary_search(target, L // 2, L)

def binary_search(target, lowbound, highbound):
  if highbound == lowbound + 1:
    return highbound
  print("Starting binary search from %s to %s" % (lowbound, highbound))
  mid = (lowbound + highbound) // 2
  print("Searching %s" % (mid))
  midprob = get_computed_prob(mid)
  if midprob >= target:
    return binary_search(target, lowbound, mid)
  else:
    return binary_search(target, mid, highbound)

# Probability that this extreme result would occur if probability is 0.5
def accuracy_test(target, length):
  prob, successes, trialsize = trial_cache[length]
  alt = "greater" if prob > target else "less"
  return stats.binom_test(successes, n=trialsize, p=target, alternative=alt)

def monte_carlo_trial(length, trialsize, force=False):
  cached_trials = trial_cache[length][2] if length in trial_cache else 0
  if trialsize > cached_trials and not force:
    print("Computing monte carlo trial with length %s (%s trials)" % (length, trialsize))
    res = brute_probability(length, trialsize - cached_trials)
    if length in trial_cache:
      trial_cache[length][1] += res[1]
      trial_cache[length][2] += res[2]
      trial_cache[length][0] = trial_cache[length][1] / trial_cache[length][2]
    else:
      trial_cache[length] = list(res)
    save_trial_cache()
  elif trialsize > cached_trials:
    print("Computing monte carlo trial with length %s (%s trials)" % (length, trialsize))
    res = brute_probability(length, trialsize)
    trial_cache[length] = list(res)
    save_trial_cache()
  return trial_cache[length]

def monte_carlo_progression(target, certainty=0.001):
  lowbound, highbound = monte_carlo_bounds(target, 1000, certainty)

  res = monte_carlo_binary_search(target, lowbound, highbound, certainty)
  return res

def monte_carlo_verified_trial(target, length, certainty, startsize=100):
  monte_carlo_trial(length, startsize)
  while (testres := accuracy_test(target, length)) > certainty:
    print("Only %s%% (requires %s%%) sure that probability for length %s is %s than %s" % ((1 - testres) * 100, (1 - certainty) * 100, length, "greater" if trial_cache[length][0] > target else "less", target))
    monte_carlo_trial(length, trial_cache[length][2] * 2)
  return trial_cache[length]

def monte_carlo_binary_search(target, lowbound, highbound, certainty):
  if lowbound + 1 == highbound:
    return highbound
  mid = (lowbound + highbound) // 2
  midprob = monte_carlo_verified_trial(target, mid, certainty)[0]
  if midprob >= target:
    return monte_carlo_binary_search(target, lowbound, mid, certainty)
  else:
    return monte_carlo_binary_search(target, mid, highbound, certainty)
    
def monte_carlo_bounds(target, trialsize, certainty):
  low = [x for x in [1]][0]
  while True:
    res = monte_carlo_verified_trial(target, low, certainty)
    if res[0] >= target:
      return low // 2, low
    low *= 2

if __name__ == "__main__":
  print(monte_carlo_progression(0.99, 0.1))
  print(monte_carlo_progression(0.99, 0.05))
  print(monte_carlo_progression(0.99, 0.01))