import subprocess
import random
import argparse
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

N = 5

def get_args_parser():
  parser = argparse.ArgumentParser()
  parser.add_argument('--prof_mult', '-profm', action='store_true',
                      help='Enable profiling of multiplication of two square\
                      matrices')
  parser.add_argument('--dim_mult', '-dimm', type=int, default=400,
                      help='The dimension of the matrices, only if specified\
                      -profm')
  parser.add_argument('--prof_vector', '-profv', action='store_true',
                      help='Enable profiling of multiplication of square\
                      matrix and vector')
  parser.add_argument('--dim_vector', '-dimv', type=int, default=1500,
                      help='The dimension of the matrix and vector, only if\
                      specified -profv')
  parser.add_argument('--prof_lu', '-proflu', action='store_true',
                      help='Enable profiling of LU decomposition of square\
                      matrix')
  parser.add_argument('--dim_lu', '-dimlu', type=int, default=600,
                      help='The dimension of the matrix, only if specified\
                      -proflu')
  parser.add_argument('--prof_inv', '-profi', action='store_true',
                      help='Enable profiling of inversion of square matrix')
  parser.add_argument('--dim_inv', '-dimi', type=int, default=600,
                      help='The dimension of the matrix, only if specified\
                      -profi')
  return parser

def prof(dim, axe, inp, prog1, prog2, title):
  times = []
  times_par = []
  thNums = [pow(2, i) for i in range(N)]
  for th in thNums:
    process = subprocess.run(prog1, input = inp,
                             capture_output=True, text=True)
    times.append(float(process.stdout))
    process = subprocess.run(prog2, input = str(th)+' '+inp,
                             capture_output=True, text=True)
    times_par.append(float(process.stdout))
  print(title + f" for {' '.join(str(x) for x in thNums)} threads")
  print("  Simple Times: ", times)
  print("Parallel Times: ", times_par)
  print()
  axe.plot(thNums, times_par)
  axe.plot(thNums, times)
  axe.set_title(title)
  axe.set_xlabel('Threads')
  axe.set_ylabel('Time, s')
  axe.legend(['Simple', 'Parallel'])

def gen_input(dim, k):
  data = [str(random.randint(-1000, 1000)) for _ in range(k)]
  return str(dim)+ ' ' + ' '.join(data)

def main():
  parser = get_args_parser()
  args = parser.parse_args()
  print(args)
  num = sum([args.prof_mult, args.prof_vector, args.prof_lu, args.prof_inv])
  if num == 0:
    return
  fig, axis = plt.subplots(num)
  fig.set_figwidth(8)
  fig.set_figheight(5*num)
  if num == 1:
    axis = [axis]
  cur = 0;
  if args.prof_mult == True:
    title1 = f"The multiplication time of two square matrices of dimension {args.dim_mult}"
    inp = gen_input(args.dim_mult, 2 * pow(args.dim_mult, 2)) 
    prof(args.dim_mult, axis[cur], inp, 'programs/mult.out',
         'programs/pmult.out', title1)
    cur+=1
  if args.prof_vector == True:
    title2 = f"The multiplication time of square matrix and vector of dimension {args.dim_vector}"
    inp = gen_input(args.dim_vector, pow(args.dim_vector, 2) + args.dim_vector)
    prof(args.dim_vector, axis[cur], inp, 'programs/mult_vector.out',
         'programs/pmult_vector.out', title2)
    cur+=1
  if args.prof_lu == True:
    title3 = f"The LU-decomposition time of square matrix of dimension {args.dim_lu}"
    inp = gen_input(args.dim_lu, pow(args.dim_lu, 2))
    prof(args.dim_lu, axis[cur], inp, 'programs/lu.out', 'programs/plu.out',title3)
    cur+=1
  if args.prof_inv == True:
    title4 = f"The inversion time of square matrix of dimension {args.dim_inv}"
    inp = gen_input(args.dim_inv, pow(args.dim_inv, 2))
    prof(args.dim_inv, axis[cur], inp, 'programs/inverse.out',
         'programs/pinverse.out', title4)
    cur+=1

  plt.savefig('results.pdf')

main()

