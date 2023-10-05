%chk=ODENICOR_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
C  1.2126270483829236  2.7323717085757027  -0.12353949141234724 
H  0.39094953572656643  3.080810304905273  0.201171379644906 
C  2.109193632335051  3.5928376651443754  -0.6955863060410117 
H  1.9298162237164163  4.525157359979608  -0.7196797055438496 
C  3.29205874344173  3.0847414658772694  -1.2439975761813498 
H  3.9098454278406294  3.6578290193568153  -1.681387972310737 
C  3.5379725660940533  1.7324357826847343  -1.1353381953912536 
H  4.334906202305811  1.3627089552913512  -1.4972063241065814 
C  2.6161460901765388  0.9138896719817058  -0.49482566393660776 
C  0.9270415123895959  -2.200077354121198  -1.5496832588958007 
C  1.4232382027165835  -1.7173395218853256  -2.771412390019698 
H  2.0335494242490313  -0.9891580502051851  -2.801244540014739 
C  0.9842779216435763  -2.3517499108559887  -3.9554878268084495 
H  1.327740543016767  -2.0668313814933033  -4.794383916527113 
C  0.05762031079036434  -3.3870932211618756  -3.913934462893285 
H  -0.23388734750461926  -3.799207400888589  -4.7187763840328145 
C  -0.43782735138111195  -3.8132900397135776  -2.6984824910304295 
H  -1.073042416303842  -4.519880847881061  -2.668936406544631 
C  -0.01329314843309004  -3.2131307300137824  -1.5098582259626239 
H  -0.3679251021544647  -3.501002284821272  -0.6766897584048976 
C  2.0775400894228673  -2.6963860262396766  1.0944316527947975 
C  2.321953263545985  -2.362324163109622  2.4184985400467816 
H  2.0970077732781762  -1.497366976550104  2.7375402486480525 
C  2.8928840507974494  -3.2863594068083444  3.277585214702532 
H  3.063895820694401  -3.055373144489816  4.183448808836238 
C  3.210436300844097  -4.534976924335289  2.8112158831474123 
H  3.6097217209136048  -5.165761144344722  3.4003417365157933 
C  2.9618562959232504  -4.896074709818602  1.504021352684327 
H  3.1691677885620346  -5.773528330008169  1.2074500974418454 
C  2.405287728763639  -3.9717384729814844  0.6182420566293192 
H  2.2523994230123887  -4.203053805627422  -0.29079829057614287 
C  6.628676514098046  -1.8734408510458025  0.43452546050658136 
H  7.496136656436677  -2.180373740047693  0.6714073574570925 
C  5.9222294451381075  -2.4782259485377347  -0.5674980594816512 
H  6.290770162514788  -3.2184602768026886  -1.034980528557881 
C  4.654030568702777  -1.9968994346421964  -0.8962164447742839 
H  4.146454697505639  -2.3791545306602426  -1.6035060082111485 
C  4.167244891617421  -0.9385339779566659  -0.14580540762238003 
C  6.058046815734011  -0.8121443981645263  1.0892260000256795 
H  6.570211722313627  -0.38177474919952137  1.7632382356105125 
N  2.8145186361983363  -0.46767868153931735  -0.3676389211760541 
N  1.445987071173182  1.410397599969597  -2.6747572984751435e-15 
N  4.817668011686986  -0.3391346383532634  0.8408704492637569 
P  1.4459870711731833  -1.410397599969597  0.0 
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


