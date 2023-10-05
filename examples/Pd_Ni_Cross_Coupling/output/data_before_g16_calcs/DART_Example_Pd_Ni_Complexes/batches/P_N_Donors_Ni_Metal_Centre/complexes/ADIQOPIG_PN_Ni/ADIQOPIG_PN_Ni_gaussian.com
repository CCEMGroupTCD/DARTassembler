%chk=ADIQOPIG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4805949201250146  1.3740228100362826  -1.3824820674711525e-16 
N  1.4805949201250148  -1.3740228100362821  2.220446049250313e-16 
N  1.5806511487145187  -2.681785422763614  -0.32149107277345373 
C  2.825714322979678  -3.1122010284436863  -0.11188508889930293 
C  3.5654805366402744  -2.069042807538932  0.34050789876199516 
C  2.6987497337069577  -0.9864159296338109  0.39836477807637655 
C  2.8606145271142003  0.4517702457091002  0.7930911056854744 
C  1.2164364156537886  2.893107045577382  0.931772241401329 
C  0.859777045300744  4.055754749335531  0.26629345099153345 
C  0.6224119149663202  5.213019038097965  0.9698261445150909 
C  0.7098414156454098  5.217650646069476  2.3294467977190045 
C  1.0418230939127555  4.076431261975994  2.997610477264754 
C  1.2967923332967515  2.902030052466676  2.3031039286483788 
C  2.076038615699214  1.8812728133907357  -1.6221224890549755 
C  1.380702420838586  1.506321112271308  -2.767327168553125 
C  1.861535440035374  1.8604686755036473  -4.013124604377543 
C  3.01209172758289  2.6033209466191805  -4.128362933432655 
C  3.69719785689237  2.9762900605645792  -2.9928461384774145 
C  3.240093600700482  2.6196060620019512  -1.756187756327897 
H  0.9354790548344218  -3.1666082729624474  -0.6177749116905615 
H  3.1331372821594563  -3.97830586040683  -0.2514438013441017 
H  4.467583259659923  -2.0778575811449684  0.5649600025320698 
H  3.715201395223412  0.7934436091866712  0.48592219849398044 
H  2.8185331842788175  0.5439770664944519  1.757671845067024 
H  0.7815925235935491  4.054254439453523  -0.6601598611811295 
H  0.40101199977169744  5.993904735309695  0.5162143882278732 
H  0.5436328640109944  6.000860744285809  2.8025017833343675 
H  1.0979777652220302  4.084586082928132  3.925892224249683 
H  1.5200929233526124  2.1256754117661827  2.7631745518613 
H  0.5925170108721649  1.0171118028696324  -2.6944166678634813 
H  1.4050894635018851  1.593575897362419  -4.777781687966348 
H  3.3253817885155037  2.8530817529852253  -4.967317030701563 
H  4.477793382624617  3.4743257845726343  -3.069360534630041 
H  3.7151657332253083  2.874348026577312  -0.997820668182184 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz


