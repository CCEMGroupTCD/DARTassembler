%chk=IVOSAJAG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
N  1.508130644208252  -1.489174925923748  -5.944053462456015e-15 
C  1.759716301340316  -2.094631158300688  1.3391631072432497 
H  2.5044631827154924  -2.727553767721698  1.2762769945818753 
H  0.954031020408401  -2.565756643820662  1.6400146071285624 
H  1.9842187101912256  -1.3884273526277928  1.9797529239538096 
C  1.1665849008063283  -2.580996156388955  -0.959849678788031 
H  0.3986216011921544  -3.0856741273357216  -0.6197805445126919 
H  1.934208558514865  -3.1817230048114666  -1.0591250217012376 
H  0.9413436912365897  -2.1916151173278045  -1.8305667476061642 
C  2.755352506174182  -0.8342063772188512  -0.4975539038553403 
H  2.671950336874172  -0.6727662295136576  -1.4696682936598515 
H  3.5210482511237258  -1.4455442020450902  -0.3554199877660048 
C  3.028329420075355  0.4791706115769887  0.20622062664329943 
H  3.805886280037387  0.9366694910839654  -0.2006659973392136 
H  3.2165386787559953  0.327773756798396  1.1661926821272424 
P  1.5081306442082536  1.4891749259237488  -2.2400469406374495e-16 
C  1.6969773623943927  2.245868394490583  -1.642167020025273 
C  0.7809241661211247  1.941705577809373  -2.650541123493907 
H  0.02581835755965156  1.4002269469495299  -2.454318596247699 
C  0.9656080715016244  2.426653442179323  -3.939297762249438 
H  0.3404812854659258  2.210640874145747  -4.621658848355797 
C  2.0468105025899175  3.213355406875166  -4.228402953749856 
H  2.17214046294983  3.539195744066585  -5.1115314697061605 
C  2.967573014146871  3.5399690581883525  -3.2265374886861933 
H  3.710340847400473  4.098355709430275  -3.424844053374813 
C  2.7891946495486364  3.047048228525131  -1.9505707670482688 
H  3.4224504706990437  3.2578890737033674  -1.2737803246041575 
C  1.6443860685817238  2.8145479131003155  1.2356531447746582 
C  1.722945161932429  4.177726532791866  0.9260242222435747 
H  1.6487024383224345  4.468009354826255  0.023858455201433808 
C  1.9096848039742913  5.1095498101272865  1.950491359786785 
H  1.9765371051625422  6.035055428041894  1.74464371944487 
C  1.995863252536536  4.689383455518629  3.252123992083013 
H  2.145398784664198  5.328137092476324  3.9391528240140903 
C  1.8674793274882189  3.3462204159185247  3.5862397008448106 
H  1.89708978151248  3.0717636325782656  4.494306909310673 
C  1.6987411257728708  2.4118827395458795  2.576736384459729 
H  1.6189315736484535  1.491405718251882  2.7962545355530573 
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
