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

# FREQUENCIES = {
#   "A": 0.4,
#   "B": 0.3,
#   "C": 0.2,
#   "D": 0.05,
#   "E": 0.025
# }

# PREFERENCES = ["A", "B"]

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
    # print(prefsremaining, (length + prefsremaining * (prefsremaining + 3) // 2) // (prefsremaining + 2) + 1, math.ceil((length - 1) / (FREQSIZE - index + 1)) + 1)
    mintodistributeleftovers = math.floor((length + prefsremaining * (FREQSIZE + (1 - prefsremaining) / 2)) / (FREQSIZE + 1)) + 1
    # lowbound = max(1 + prefsremaining, mintodistributeleftovers, math.ceil((length - 1) / (FREQSIZE - index + 1)) + 1)
    lowbound = max(1 + prefsremaining, mintodistributeleftovers)
    # At minimum, remaining n prefs need 1 + 2 + 3 + ... + n (assuming strict inequality)
    highbound = length - (prefsremaining) * (prefsremaining + 1) // 2
    if counts:
      highbound = min(highbound, counts[-1] - 1)
    # print(counts, index, prefsremaining, length, lowbound, highbound)
  else:
    lettersafter = FREQSIZE - index
    lowestpref = counts[PREFSIZE - 1]
    lowbound = max(0, length - lettersafter * (lowestpref - 1))
    highbound = min(length, lowestpref - 1)

    cache_key = (lowestpref, length, index)
    if cache_key in sum_cache:
      return sum_cache[cache_key]
    # cache_key = (index, lowestpref, length)
    # subcache_key = (lowbound, highbound)
    # if cache_key in sum_cache and subcache_key in sum_cache[cache_key]:
    #   return sum_cache[cache_key][subcache_key]
  # print(counts, index, lowbound, highbound)

  # prob = D(ORDERED_FREQUENCIES[index])
  tot = 0
  # if index >= PREFSIZE and (cache_key in sum_cache) and (closest := closest_cache(sum_cache[cache_key], subcache_key)):
  #   tot = sum_cache[cache_key][closest]
  #   lowmidbound = closest[0] - 1
  #   highmidbound = closest[1] + 1
  #   r = it.chain(range(lowbound, lowmidbound + 1), range(highmidbound, highbound + 1))
  #   # print(closest, (lowbound, highbound), [i for i in it.chain(range(lowbound, lowmidbound + 1), range(highmidbound, highbound + 1))])
  #   print(closest, (lowbound, highbound))
  # else:
  #   r = range(lowbound, highbound + 1)
  if index == 0:
    bar = alive_bar(highbound - lowbound + 1)
  else:
    bar = nullcontext()
  with bar as bar:
    for c in range(lowbound, highbound + 1):
      # print(c)
      nextiter = recursive_sum(length - c, index + 1, counts + [c]) if index + 1 < FREQSIZE else 1
      toadd = nextiter * exponent(index, c) / factorial(c)
      if index + 1 == FREQSIZE:
        # print("Reached end:", counts + [c, length - c])
        toadd *= exponent(-1, length - c) / factorial(length - c)
      tot += toadd
      if index == 0:
        bar()
  if index >= PREFSIZE:
    sum_cache[cache_key] = tot
    # if cache_key not in sum_cache:
    #   sum_cache[cache_key] = {}
    # sum_cache[cache_key][subcache_key] = tot
  # print("Returning")
  return tot

# factorial_cache = { 0: 1 }
# def factorial(x):
#   if x not in factorial_cache:
#     f = 1
#     for i in range(x, 1, -1):
#       if i in factorial_cache:
#         f *= factorial_cache[i]
#         break
#       else:
#         f *= i
#     factorial_cache[x] = f
#   return factorial_cache[x]

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
  
# def probability(n, trials):
#   success = 0
#   tot = 0
#   for i in range(trials):
#     tot += 1
#     string = [random_letter() for i in range(n)]
#     counts = {}
#     E = len([c for c in string if c == "E"])
#     T = len([c for c in string if c == "T"])
#     if T >= E:
#       continue
#     works = True
#     for l in FREQUENCIES.keys():
#       if l != "T" and l != "E":
#         count = len([c for c in string if c == l])
#         if count >= T or count >= E:
#           works = False
#           break
#     if works:
#       success += 1
#   return success / tot

# def min_len_for(p, trials):
#   if p < 0 or p > 1: return False
#   min = 1
#   max = 1000
#   while True:
#     n = (min + max) // 2
#     prob = probability(n, trials)
#     print(n, prob)
#     if prob >= p:
#       max = n
#     if prob < p:
#       min = n
#       probmax = probability(max, trials)
#       if probmax < p:
#         min = max
#         max *= 2
#     if min + 1 == max:
#       if prob < p:
#         max *= 2
#       else:
#         return min
#     n += 1

# def min_len_for(p, trials):
#   n = 1
#   while True:
#     prob = probability(n, trials)
#     print(n, prob)
#     if prob >= p:
#       return n
#     n += 1

# L = 100
# print(brute_probability(L, 10000))
# print(get_computed_prob(L))

print(get_computed_prob(100))
# print(brute_probability(618, 1000000))
# print(brute_probability(619, 1000000))

# print(compute_prob(100))
# print(smart_search(0.5))

# for i in range(100, 100 * 100, 100):
#   print("%s - %s" % (i, get_computed_prob(i)))

# for i in range(26, 100):
#   brute = brute_probability(i, 10000)
#   comp = float(compute_prob(i))
#   err = abs((brute - comp) / (brute) * 100)
#   print("Length: %s Brute: %s Calculated: %s Error: %s%%" % (i, brute, round(comp, 4), round(err, 2)))

# l = min_len_for(0.5, 100)
# print(probability(l, 10000))

# 100 -> 0.105381
# 300 -> 0.300386