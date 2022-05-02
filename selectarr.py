probs = [
  1202,  910,  812,  149,  271,
  432,  230,  203,  592,  731,
  10,   69,  398,  261,  695,
  768,  182,   11,  602,  628,
  288,  111,  209,   17,  211,    8 ]

arr = []
for i in range(len(probs)):
  arr += [str(i + 1).zfill(2)] * probs[i]

# print(arr)
lines = ["    integer, dimension(10000) :: selectarr = (/ &"]
linebase = "      "
linesubstrstart = 0
linesubstrlen = 1
lastfullline = ""
while True:
  fullline = linebase + ", ".join(arr[linesubstrstart:linesubstrstart+linesubstrlen]) + ", &"
  if len(fullline) > 132:
    lines.append(lastfullline)
    linesubstrstart += linesubstrlen - 1
    linesubstrlen = 1
  if linesubstrlen + linesubstrstart > 10000:
    lines.append(linebase + ", ".join(arr[linesubstrstart:linesubstrstart+linesubstrlen]) + " /)")
    break
  lastfullline = fullline
  linesubstrlen += 1

with open("selectarr.txt", "w+") as f:
  f.write("\n".join(lines))