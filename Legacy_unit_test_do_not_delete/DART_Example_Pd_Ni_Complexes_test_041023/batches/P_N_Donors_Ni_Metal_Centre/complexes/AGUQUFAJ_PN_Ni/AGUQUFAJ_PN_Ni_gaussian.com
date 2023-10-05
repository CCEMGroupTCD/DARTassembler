%chk=AGUQUFAJ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
C  1.2867011104523371  -2.7323717085757027  0.2105870102649476 
H  0.4047796384044504  -3.080810304905273  0.2635549414005571 
C  2.341023817517669  -3.5928376651443754  0.35013208166777193 
H  2.188634985885459  -4.525157359979608  0.4477762858793377 
C  3.6448322964333437  -3.0847414658772694  0.34726108916499554 
H  4.3895863481458965  -3.6578290193568153  0.4825834764648909 
C  3.821784470211762  -1.7324357826847343  0.1448545739530151 
H  4.6969837099989435  -1.3627089552913512  0.13601976882958872 
C  2.7156336641246686  -0.9138896719817058  -0.046066218120072666 
C  1.6305871153227387  2.200077354121198  1.6238058749687725 
C  2.596619081381458  1.7173395218853258  2.5213667174084016 
H  3.1623565053470575  0.9891580502051854  2.290475159688066 
C  2.699197863310717  2.351749910855989  3.780016137169152 
H  3.365013518835232  2.0668313814933037  4.395160619778131 
C  1.841799644250368  3.387093221161876  4.133978428602845 
H  1.9177448772448389  3.7992074008885894  4.9866093888978655 
C  0.8790993703772694  3.813290039713578  3.2417900715048913 
H  0.290912295728176  4.519880847881061  3.483465711628119 
C  0.7615237036066861  3.213130730013783  1.9851147373138571 
H  0.08800579292057464  3.501002284821272  1.3798816070923334 
C  1.5558416868412328  2.6963860262396766  -1.25879776807709 
C  1.2177804034845066  2.3623241631096215  -2.5621033693453805 
H  0.8790777016841622  1.4973669765501036  -2.7561872821567985 
C  1.3721537047018617  3.286359406808344  -3.581986089207771 
H  1.1443085059074187  3.0553731444898156  -4.475250015495207 
C  1.8570499779411445  4.534976924335289  -3.293515312310851 
H  1.9699501193032878  5.165761144344721  -3.996189970931407 
C  2.1842042640010324  4.896074709818602  -2.003740280597297 
H  2.497428712494589  5.773528330008169  -1.8225690652176751 
C  2.054128343949114  3.9717384729814844  -0.9657355665942946 
H  2.29974153328805  4.203053805627422  -0.07725183117441049 
C  5.959460476287255  1.8734408510458023  -2.5841130120335825 
H  6.645535742705613  2.1803737400476924  -3.165405417332337 
C  5.742675434800369  2.4782259485377347  -1.3774138607787119 
H  6.274253385274174  3.2184602768026886  -1.1094828961102083 
C  4.732219310170481  1.9968994346421964  -0.5435298237452771 
H  4.57111283143963  2.3791545306602426  0.3120030476741615 
C  3.9739042522979418  0.9385339779566658  -1.0179086735042135 
C  5.165605933173296  0.812144398164526  -2.936314677956609 
H  5.344935096858475  0.3817747491995209  -3.7636274581683886 
N  2.8416688071765828  0.4676786815393173  -0.24517241406859888 
N  1.4459870711731833  -1.410397599969597  3.1387312607223207e-15 
N  4.146400555346927  0.33913463835326313  -2.18702137428094 
P  1.4459870711731833  1.410397599969597  -1.727238906327876e-16 
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


