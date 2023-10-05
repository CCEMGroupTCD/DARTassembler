%chk=NIKOKIJO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3122583453725871  -1.5356034758361279  5.41489398294622e-16 
N  1.3122583453725873  1.5356034758361283  3.670543241683518e-16 
N  2.9334068644896254  -1.0433422300816488  -0.24712598991800624 
C  2.590499086576874  1.4113331951270451  0.14102143403769882 
H  3.052054486375652  2.198640686301946  0.40821639517371533 
C  3.402103049004936  0.2590308606741194  -0.04909147500694613 
C  4.777890311546627  0.2173359010002769  -0.10048007846961177 
H  5.360399684395281  0.9589777590195757  0.005089400913200031 
C  5.16918351349805  -1.1070075489861706  -0.3369435990908908 
H  6.061456889710804  -1.420635179802626  -0.4351432492448116 
C  4.04450874637635  -1.8573647851682553  -0.40024799316399823 
H  4.020759399102035  -2.797343456559356  -0.5303538186825572 
C  0.774160802028074  2.8382373793644007  0.3295138346493215 
C  1.0299747581131427  3.9240585670848267  -0.5243941394143484 
C  0.5149315720240532  5.158885245871501  -0.14918066708176025 
H  0.6835489255853641  5.915033947845336  -0.6979961721299249 
C  -0.2304053549012992  5.3122952919589075  0.992142894893974 
H  -0.5647467438940923  6.169562412485066  1.2308185318390958 
C  -0.4936687317110331  4.216365371584336  1.7969749109436264 
H  -1.0150214975547498  4.329025181048367  2.5817778931491366 
C  -0.011684250135256757  2.9589955489817426  1.4785189819228748 
C  1.7792027961052816  3.788617271025794  -1.8351156744508144 
H  1.9921320540798726  2.8198248168475017  -1.9720934530340108 
C  3.1014855858304857  4.561458739331151  -1.7926386351920025 
H  2.9183035795872527  5.515940018557176  -1.6580797997490697 
H  3.5792789999357053  4.436504818596163  -2.6389593197617813 
H  3.6518569164040118  4.226217267101676  -1.0544592845053464 
C  0.9219617456380466  4.253087963638444  -3.01221100070855 
H  0.08205679443313163  3.7501387178512338  -3.022177901106999 
H  1.4066239336226123  4.097469567509273  -3.849531135617672 
H  0.729247994583311  5.208017193791143  -2.919287669115291 
C  -0.31449112079250297  1.7561796631982787  2.3572979552658424 
H  -0.29327419826240386  0.9475068452861904  1.7717212773077182 
C  0.7559201382050066  1.5758244677600817  3.411960171375763 
H  1.629041462689259  1.474968512797053  2.9781124572185873 
H  0.5608703169817071  0.7761243209944279  3.9422630777450887 
H  0.7719037506018673  2.361502264134613  3.9993142591446325 
C  -1.697619908412948  1.8171372290497305  2.9962734493351544 
H  -1.7399882709787153  2.5791811247649528  3.612457482604991 
H  -1.8664054824142924  0.988457732666644  3.491196324421078 
H  -2.3759260916513796  1.9238808758174932  2.296657257659144 
C  1.1848277458210192  -1.9123113917159247  1.7633580119792704 
C  2.3119010902917436  -2.0940947618323578  2.555120364287548 
H  3.1798169655411055  -2.0327151632034344  2.173655361486115 
C  2.1590401074884777  -2.3660001310049905  3.90515258247775 
H  2.927050870906325  -2.4920890224713066  4.4496257762377125 
C  0.906550464263571  -2.4569009020676473  4.464145519315618 
H  0.8135483928232488  -2.645167983738965  5.390949905724522 
C  -0.21410300676483907  -2.2748605971329448  3.6816750564537197 
H  -1.0789721898847282  -2.3368307225575333  4.071541186868281 
C  -0.08227186670675724  -2.0004419880121547  2.3320779693646645 
H  -0.8549856851692135  -1.8738656041613917  1.7939713778714352 
C  1.2736872102656374  -3.0863006821597954  -0.9100033032008982 
C  1.044411141542033  -4.286740143913605  -0.26557053481083237 
H  0.9039366736297314  -4.305790878989644  0.6722229197183821 
C  1.0202413571715452  -5.469295757113967  -0.9991902545647463 
H  0.8610933778245418  -6.29720268100769  -0.5625591447722804 
C  1.2278566770044712  -5.434889107981295  -2.3569296483609943 
H  1.2090883788895384  -6.243492217142599  -2.856918551078379 
C  1.464494597955372  -4.231721139118763  -3.0056342290996025 
H  1.6172778143296245  -4.2185197493596265  -3.942756959643886 
C  1.476694882429896  -3.0615202941242208  -2.2896786191398726 
H  1.6232617102081475  -2.2349794422881653  -2.733678478645648 
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

