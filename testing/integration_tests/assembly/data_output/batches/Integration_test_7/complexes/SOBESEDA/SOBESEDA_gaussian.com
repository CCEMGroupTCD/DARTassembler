%chk=SOBESEDA_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program

-1 0
Fe  0.0  0.0  0.0 
N  1.51032  1.10265  0.0 
N  1.51032  -1.10265  0.0 
C  2.40643  -0.00805  -0.24581 
C  3.88372  0.0349  0.30285 
C  4.58507  -1.32553  0.15057 
H  5.66776  -1.24134  0.38248 
H  4.48634  -1.71688  -0.87549 
H  4.18622  -2.06401  0.8546 
C  4.66288  1.0212  -0.58932 
H  4.32492  2.06104  -0.44126 
H  4.53376  0.74389  -1.65776 
H  5.74799  1.01066  -0.36325 
C  4.00861  0.48265  1.76351 
H  3.42449  -0.15871  2.4485 
H  3.67578  1.52187  1.85957 
H  5.07363  0.46553  2.08322 
C  1.7955  2.46008  -0.32043 
C  1.88797  3.42129  0.72753 
C  2.35153  4.71827  0.42445 
H  2.43325  5.43737  1.21711 
C  2.75744  5.06819  -0.87591 
C  2.59397  4.11829  -1.89251 
H  2.87006  4.40856  -2.89137 
C  2.12865  2.82744  -1.65062 
C  1.52671  3.06001  2.17503 
H  1.53324  1.95645  2.30039 
C  0.10562  3.51529  2.47114 
H  0.00606  4.61154  2.34835 
H  -0.19356  3.20883  3.49445 
H  -0.58366  3.03176  1.77268 
C  2.47932  3.65072  3.23682 
H  3.53537  3.36997  3.05672 
H  2.18699  3.28073  4.24287 
H  2.41309  4.75738  3.27406 
C  2.01515  1.84504  -2.81024 
H  2.27433  0.8474  -2.46224 
C  3.01352  2.08084  -3.9537 
H  2.97451  1.22766  -4.66235 
H  4.04609  2.1457  -3.54568 
H  2.78279  3.00302  -4.52277 
C  0.57775  1.83113  -3.32007 
H  -0.12071  1.60265  -2.48908 
H  0.45652  1.05678  -4.10409 
H  0.31811  2.82772  -3.73866 
C  3.36136  6.43939  -1.27437 
C  4.83192  6.27943  -1.75328 
C  5.56493  7.42221  -2.13933 
H  5.1006  8.40378  -2.09912 
C  6.88601  7.32663  -2.58023 
H  7.42541  8.21832  -2.87406 
C  7.51252  6.08298  -2.64611 
H  8.53401  6.00521  -2.99328 
C  6.82495  4.93764  -2.25635 
H  7.31068  3.96985  -2.29723 
C  5.51125  5.03155  -1.8089 
H  5.04663  4.10438  -1.5064 
C  2.43949  6.92653  -2.42549 
C  2.90721  7.15991  -3.74363 
H  3.94463  7.0677  -4.0178 
C  2.02309  7.48345  -4.77146 
H  2.40202  7.6419  -5.77215 
C  0.65652  7.57642  -4.51913 
H  -0.02542  7.80639  -5.32502 
C  0.16981  7.35305  -3.2322 
H  -0.89154  7.41065  -3.03786 
C  1.04899  7.03281  -2.19602 
H  0.64207  6.83582  -1.21185 
C  3.44049  7.41705  -0.08261 
C  2.70463  8.62558  -0.01632 
H  2.02362  8.92462  -0.79645 
C  2.84982  9.49159  1.06723 
H  2.28132  10.41405  1.10209 
C  3.73397  9.18951  2.10237 
H  3.85004  9.86999  2.93638 
C  4.47548  8.00552  2.05804 
H  5.16972  7.77522  2.85512 
C  4.33501  7.12884  0.97867 
H  4.92258  6.22726  0.97175 
C  1.8003  -2.48857  -0.1183 
C  2.01333  -3.0391  -1.40149 
C  2.4128  -4.37378  -1.51271 
H  2.58539  -4.73792  -2.50614 
C  2.60362  -5.20248  -0.38712 
C  2.37644  -4.63777  0.87649 
H  2.54223  -5.25775  1.74635 
C  1.98593  -3.29707  1.03603 
C  1.81234  -2.21803  -2.67339 
H  1.52494  -1.18798  -2.41912 
C  3.07895  -2.1476  -3.54432 
H  3.92926  -1.71138  -2.98833 
H  2.89265  -1.49916  -4.426 
H  3.37537  -3.15106  -3.9147 
C  0.64345  -2.78065  -3.48011 
H  0.45576  -2.13844  -4.36622 
H  -0.27165  -2.78648  -2.84984 
H  0.85419  -3.81453  -3.82354 
C  1.80903  -2.74355  2.44416 
H  1.70123  -1.64169  2.40323 
C  0.51771  -3.29631  3.03489 
H  -0.33388  -2.99979  2.38669 
H  0.35289  -2.87858  4.04949 
H  0.55566  -4.40461  3.09455 
C  3.00096  -3.06801  3.37615 
H  3.08407  -4.15798  3.57229 
H  2.85828  -2.5731  4.35952 
H  3.96627  -2.72119  2.95893 
C  3.11217  -6.66991  -0.50373 
C  3.45231  -6.99074  -1.98465 
C  2.40074  -7.17943  -2.91371 
H  1.36833  -7.09419  -2.59599 
C  2.66631  -7.48388  -4.24774 
H  1.84654  -7.63152  -4.93728 
C  3.98284  -7.59658  -4.69106 
H  4.1871  -7.82813  -5.72711 
C  5.03714  -7.41041  -3.79906 
H  6.0601  -7.49232  -4.1427 
C  4.78126  -7.11926  -2.45945 
H  5.63864  -6.98801  -1.81881 
C  4.39744  -6.81326  0.37242 
C  5.34299  -5.75715  0.42998 
H  5.17566  -4.83265  -0.10692 
C  6.51897  -5.8788  1.17623 
H  7.22571  -5.06063  1.20824 
C  6.79249  -7.05623  1.86953 
H  7.70987  -7.15677  2.43416 
C  5.88432  -8.10876  1.83291 
H  6.0984  -9.02769  2.36437 
C  4.69974  -7.99089  1.10604 
H  4.05625  -8.84963  1.12584 
C  2.02308  -7.69197  -0.06216 
C  2.16967  -9.05641  -0.4138 
H  3.03497  -9.38817  -0.97442 
C  1.2106  -10.00368  -0.06255 
H  1.34754  -11.03678  -0.35061 
C  0.08089  -9.62066  0.65471 
H  -0.65953  -10.35752  0.9319 
C  -0.09852  -8.28649  1.00782 
H  -0.97755  -7.98526  1.56293 
C  0.84842  -7.32925  0.64595 
H  0.62633  -6.32207  0.9445 
Cl  -1.52028  1.52028  0.0 
Cl  -1.52028  -1.52028  0.0 

-x 0
6-31+g(d)
****
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-H 0
6-31g(d,p)
****
-x 0
6-31+g(d)
****


