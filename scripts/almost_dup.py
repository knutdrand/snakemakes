from collections import deque, Counter
lines = (line.split()[1:3] for line in open(snakemake.input[0]))
intervals = ((int(start), int(end)) for start, end in lines)
max_diff = 4
counts = Counter()
nearby_starts = deque([])
t = 0
for start, end in intervals:
    if nearby_starts and start<nearby_starts[-1][0]:
        nearby_starts = deque([])
    while nearby_starts and (start-nearby_starts[0][0]>max_diff):
        nearby_starts.popleft()
    diffs = ((start-interval[0])+abs(end-interval[1]) for interval in nearby_starts)
    ds = [d for d in diffs if d<=max_diff]
    counts.update(ds)
    t+=int(bool(ds))
    nearby_starts.append((start, end))
    
print(counts)
print(t)
open(snakemake.output[0], "w").write(str(counts))
