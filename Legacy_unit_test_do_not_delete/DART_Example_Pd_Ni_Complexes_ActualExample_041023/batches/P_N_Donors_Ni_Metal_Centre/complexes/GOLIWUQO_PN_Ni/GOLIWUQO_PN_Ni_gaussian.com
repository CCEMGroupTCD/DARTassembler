%chk=GOLIWUQO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.3547405563428008  -1.4982583305291506  1.8348372687783814e-16 
N  1.3547405563428008  1.4982583305291506  -1.8348372687783814e-16 
N  2.6456602061870624  0.8711426850795698  0.12399864869568722 
C  1.4362937351413785  2.5069000708799463  -1.0960869849833976 
H  2.075439033586141  3.208508330670994  -0.8481870913389218 
H  0.5530146645180556  2.9066884103103385  -1.2410702497221653 
H  1.735915666536419  2.0696377499231713  -1.9206971115878504 
C  1.024009332917445  2.1645154302960856  1.292665105455306 
H  1.023075856233023  1.4994915591135352  2.0121460846192165 
H  0.13908637055860495  2.57965688976678  1.2295884907822672 
H  1.6949190030763046  2.853374380031066  1.486570834760441 
C  1.7226610741795176  -2.950358879992886  -0.9822944094322726 
C  1.6670176751121915  -4.017032745530859  -0.15446516573613162 
C  1.3156318446922883  -3.7027196458441844  1.2678265882811346 
C  1.1094112163161762  -2.403754133459783  1.536212392801582 
C  2.820043436580117  -0.4965946324738444  0.12095741054784476 
C  4.167243687092215  -0.7223139981593074  0.2584479301291166 
H  4.591490240878929  -1.5716931766822468  0.2924714909627583 
C  4.8046173262797875  0.5195327095061693  0.34028611638232026 
H  5.739692224186021  0.6609150170052951  0.4338543718536531 
C  3.8535126851207604  1.4909432835998744  0.26361253485189023 
H  4.002621232173453  2.427652879270428  0.30185367205438624 
C  2.010318760835426  -2.8706427572346596  -2.4451391602791297 
H  2.1510490673182368  -3.772570070078461  -2.800695348783305 
H  2.816010430382377  -2.332192929483607  -2.5901449071222284 
H  1.2515205612113929  -2.45332800432474  -2.904265318596295 
C  1.9429298954157745  -5.429629543342436  -0.5439942448812298 
H  2.167170930571956  -5.46873658103153  -1.496882159462627 
H  1.1484270492831086  -5.976719953786212  -0.3723314058986296 
H  2.694432820518533  -5.77268900413547  -0.01707582622194044 
C  1.2650943404292305  -4.795025933520227  2.2956871014667826 
H  0.9805434238087516  -4.422007747262174  3.156132206381229 
H  2.1539021639217912  -5.1962710911478025  2.3916200800448033 
H  0.6273870173775717  -5.482026021573446  2.011278119292921 
C  0.7759728962883957  -1.7275026190482703  2.821914502931788 
H  0.49538068128564117  -2.3975697460205825  3.479941876398055 
H  0.04846082146626496  -1.08845767233214  2.6762840216973145 
H  1.5650574985654389  -1.2535789100779886  3.155584382491151 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Ni 0
lanl2dz


