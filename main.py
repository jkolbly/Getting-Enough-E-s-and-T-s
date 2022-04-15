import random
import math
import decimal
from decimal import Decimal as D
from alive_progress import alive_bar
from contextlib import nullcontext
import itertools as it

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

def closest_cache(subcache, key):
  closenesses = {k:(k[0] - key[0] + key[1] - k[1]) for k in subcache.keys() if k[0] >= key[0] and k[1] <= key[1]}
  if not closenesses:
    return None
  return min(closenesses, key=closenesses.get)

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

def random_letter():
  rand = random.random()
  for l, f in FREQUENCIES.items():
    rand -= f
    if rand < 0:
      return l
  return "Z"

def test_distribution(letter, num):
  count = 0
  for i in range(num):
    if random_letter() == letter:
      count += 1
  print("Got %s for frequency of %s (expected %s)" % (count / num, letter, FREQUENCIES[letter]))

def random_index():
  rand = random.random()
  for i in range(len(ORDERED_FREQUENCIES)):
    rand -= ORDERED_FREQUENCIES[i]
    if rand <= 0:
      return i
  return len(ORDERED_FREQUENCIES)

def single_trial(length):
  counts = [0] * (len(ORDERED_FREQUENCIES) + 1)
  for i in range(length):
    counts[random_index()] += 1
  for p in range(PREFSIZE):
    c = counts.pop(0)
    if next((x for x in counts if x >= c), None) is not None:
      return False
  return True

def brute_probability(length, trials):
  success = 0
  tot = 0
  with alive_bar(trials) as bar:
    for i in range(trials):
      if single_trial(length):
        success += 1
      tot += 1
      bar()
  return success / tot

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

print("Solution: %s Probability: %s Brute forced: %s" % (sol := smart_search(0.5), get_computed_prob(sol), brute_probability(sol, 10000)))