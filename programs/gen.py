import random
import sys

N=int(sys.argv[1])

if len(sys.argv) == 3:
  print(sys.argv[2]) # number of threads

print(N)
print(' '.join([str(random.randint(1, 2)) for _ in range(N*N)]))
