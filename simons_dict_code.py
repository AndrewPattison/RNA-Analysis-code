    
sums = defaultdict(int)
counts = defaultdict(int)
positions= []
average_As = []

for x in sums.keys():
    print "average " + str(x) + ": " +  str(float(sums[x])/counts[x])
    positions.append(int(x))
    average_As.append(float(sums[x])/counts[x])



for k in newlist:
    sums[k[0]] += k[1]
    counts[k[0]] += 1
    
    
    