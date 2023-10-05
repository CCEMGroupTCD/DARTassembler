%chk=UHATEYUB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  2.7065359951413805  -1.1357138486741563  1.2590343174883658 
H  2.212190798351559  -1.1239483918710285  2.105483140771289 
C  3.2677798265729834  0.2836435684720009  1.1333602387284223 
H  3.9337716885864875  0.4007348066240435  1.828230832205958 
H  3.7292576765302368  0.3442524737243291  0.28206484663296966 
C  2.298825520338976  1.4510299079422686  1.2195694426963672 
H  2.83932315109337  2.2667223019052893  1.2544711558316166 
C  3.8074977927220495  -2.149897960394884  1.4847442370780082 
H  4.2892117897677595  -2.286107535025887  0.6671468319762714 
H  3.422338544079672  -2.9801592545998616  1.7724015792932029 
H  4.4094298399877  -1.8243224506722584  2.15950417064214 
C  2.150564907250349  1.8315727687941457  -1.3015090133738467 
H  2.6925879780918285  1.044764664653552  -1.5178754570157975 
C  3.0915176717763133  3.0311342596529363  -1.1835375567430542 
H  3.4569802826209575  3.2361670282515727  -2.047577812677576 
H  3.803723193993635  2.8197324198401574  -0.5765058802656917 
H  2.6036960455888942  3.788975379217664  -0.852724107015344 
C  1.1497050795222683  2.022570067173828  -2.419963222204852 
H  0.4675292554696616  2.637679106828455  -2.1384663937798707 
H  0.7482395090758232  1.1783703832528931  -2.639686797211965 
H  1.5965299057478428  2.374527144917631  -3.1929508159568565 
C  0.8186898604021976  -3.1264288047924547  0.6414775629812624 
C  1.0286977975327196  -4.3521644517033735  0.0029270088087209323 
H  1.4994649819248171  -4.37761041557898  -0.7982591250091322 
C  0.5524563256096807  -5.514090002432331  0.5409227011555937 
H  0.7035367346750592  -6.319908464740845  0.10194217106412534 
C  -0.14844679993547127  -5.507182386619196  1.7329578572834676 
H  -0.4807463611965377  -6.299653149340038  2.0876209394027647 
C  -0.347507513700845  -4.3280415620917765  2.3746865913822424 
H  -0.7945127323892873  -4.3185042794797805  3.1909569171357433 
C  0.11037661313061542  -3.1321803689123193  1.8271310845674664 
H  -0.060726095477488196  -2.329700731384074  2.2649715400526746 
C  1.453553517076525  1.4326527284952406  2.471748395472797 
H  0.9313315769184999  0.6271109547821427  2.49420916172755 
H  0.8684696393785732  2.19342693040884  2.475698048705639 
H  2.024585012350987  1.464290701529142  3.2438055048027543 
C  2.170679399710256  -1.9178007876969863  -1.6145034878907591 
C  3.442416489980081  -2.4340304916539646  -1.8043861779509554 
H  4.0077584441002605  -2.5481782318529866  -1.074773089856287 
C  3.876600263236897  -2.7794073603858047  -3.063965778288118 
H  4.739232460259206  -3.105210991127458  -3.182404525563094 
C  3.03995423654716  -2.6436870073246386  -4.148165982027843 
H  3.3231756346645134  -2.911261614780102  -4.992548388934526 
C  1.8038791178906761  -2.122412241893084  -3.9798163162644156 
H  1.2469626447570918  -2.003427884104873  -4.715478238993985 
C  1.3666032619454052  -1.7657464366773372  -2.719737779348064 
H  0.5117087510409616  -1.4150181591212196  -2.6167061795935216 
N  1.439777577440349  1.5553586491545919  4.0486404104429374e-17 
P  1.439777577440349  -1.5553586491545923  -5.551115123125783e-17 
H  0.9686712046314819  2.2717927015872443  -0.016231266549618373 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
****
-P 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-N 0
6-31+g(d)
****
-Br 0
6-31+g(d)
****
-C 0
6-31g(d)
****

-Pd 0
lanl2dz
