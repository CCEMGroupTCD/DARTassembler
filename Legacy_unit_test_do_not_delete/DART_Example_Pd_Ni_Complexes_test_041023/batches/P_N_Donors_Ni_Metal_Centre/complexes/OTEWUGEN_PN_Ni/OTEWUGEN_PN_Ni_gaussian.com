%chk=OTEWUGEN_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.46320914089545  1.392522534826635  0.0 
N  1.46320914089545  -1.392522534826635  0.0 
C  2.4583087780689565  -1.0064872308687556  0.7003623820532713 
C  2.6119925265816137  0.39486579660159915  1.0058577958368033 
C  3.415500022103848  0.8973306418740017  1.9627796703411886 
C  4.4082006759870715  0.21068785076180063  2.774842586126354 
C  5.279911320235941  -0.7517503909768517  2.277606461591671 
C  6.233434129945472  -1.3177095336761946  3.101404944327983 
C  6.311930589262317  -0.9627639743488303  4.416970440906773 
C  5.474210232769803  -0.02358411384997039  4.923387417272524 
C  4.536993676512995  0.5758302648221241  4.1144960900107606 
C  0.918157105875902  2.79649786815465  0.9917237783413092 
C  0.3455201707761273  2.548815046962144  2.2302779921913523 
C  -0.11623688179843272  3.598723068868602  3.005130507147664 
C  -0.04571780840847017  4.888670272648415  2.516406820644913 
C  0.5000634129876522  5.134377283804682  1.3016477074705846 
C  0.9957180082391494  4.099562611897489  0.5251286433055786 
C  2.4505674591420066  2.015500560893553  -1.3904203178992847 
C  3.8259003303446484  2.0734494399830004  -1.3444826340022047 
C  4.535796161508722  2.499341707553259  -2.45045940379429 
C  3.893972900805019  2.8753695109849566  -3.5807123325965224 
C  2.535093650981487  2.8383057804461718  -3.6402444607091633 
C  1.7970484380132556  2.39442985417475  -2.554044200309915 
H  1.4029534247544029  -2.2275207132545045  -0.19447398462126086 
H  3.0802687006969274  -1.6238392007649018  1.0108984667680174 
H  3.3242349250503325  1.8092442578056058  2.1257984702949364 
H  5.221810377750346  -1.014241827099648  1.3872851395710812 
H  6.828505303783674  -1.9436181909522574  2.755913843409692 
H  6.942247939178946  -1.3667493002093112  4.968965753281341 
H  5.534716216643016  0.21573170450516077  5.820757605656649 
H  3.980741442080687  1.232652755751415  4.4659875936043605 
H  0.2717551265275424  1.6737514917775922  2.5418954060737895 
H  -0.471556902399499  3.4377916552682475  3.8485561192064224 
H  -0.37680065137143504  5.593097608987446  3.025886848628208 
H  0.5423471522400247  6.007474587146609  0.9846981936803069 
H  1.3792514983843807  4.278009971839809  -0.3034749661885221 
H  4.27720606852068  1.8251515407730277  -0.5697929839246312 
H  5.464296238581374  2.530489935977008  -2.417689833425809 
H  4.382739547638138  3.159160659865142  -4.319491953837392 
H  2.0988988502161043  3.1121710935468165  -4.413727082405064 
H  0.8697559302376958  2.351008273203691  -2.6057525297095907 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
****
-P 0
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
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz


