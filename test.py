# python -m pip install numpy
# python -m pip install --upgrade pip
#import scipy as np
import numpy as np
import time


isLoad = True
foldername = 'D:\\_install\\_Programming\\parallel_computing\\_matrix_product_tests\\'
N_count = 20

if isLoad:
	start_time = time.time()
	print('Loading a...')
	a = np.loadtxt(foldername+'A.dat')
	print('a[0,0]=',a[0,0])
	print('Loading b...')
	b = np.loadtxt(foldername+'B.dat')
	print('b[0,0]=',b[0,0])
	print('Loading c2...')
	c2 = np.loadtxt(foldername+'C.dat')
	if not(a.shape==b.shape and a.shape==c2.shape):
		print('Error dimension of a,b and c2 are differ')
		exit()
	N = a.shape[0]
	t1 = (time.time() - start_time)
	print('Elapsed time = ',t1,'s')
	print('N=',N)
else:
	N = 2000
	a = np.random.rand(N,N)
	b = np.random.rand(N,N)
	c = np.zeros([N,N], dtype = float)


print('Matrix product in python (using numpy): c= a*b for size (',N,'x',N,')')
start_time = time.time()
for nc in range(N_count):
	c = a@b
t1 = (time.time() - start_time)/N_count
print('Elapsed time = ',t1,'s')

if isLoad:
	start_time = time.time()
	print('norm = ',np.linalg.norm(c-c2))
	t1 = (time.time() - start_time)
	print('Elapsed time = ',t1,'s')

#print('Matrix product in python: c= a*b for size (',N,',',N,')')
#start_time = time.time()
#for nc in range(N_count):
#	for i in range(N):
#		for j in range(N):
#			c3[i,j] = 0
#			for k in range(N):
#				c3[i,j] += a[i,k]*b[k,j]
#t2 = (time.time() - start_time)/N_count
#print('Elapsed time = ',t2,'s')
#print('(ratio =',t1/t2,')')

#print('norm = ',np.linalg.norm(c-c3))

