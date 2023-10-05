%chk=EPIWAMIL_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4109171476383726  -1.445480128711564  -5.45454205105122e-15 
N  2.8068489441979727  -0.46314373625926675  -0.18395081249545792 
C  4.193899727530319  -0.668220498721734  -0.2986806950677988 
H  4.620461466289198  -1.5149516060465773  -0.3535942172226277 
C  4.829828797559631  0.5332372076169904  -0.3194346771184822 
H  5.768458552201881  0.6651142125052667  -0.37722740401307986 
C  3.8384893277964625  1.561086591699715  -0.23836366299830222 
C  3.781079011847194  2.9571923630040855  -0.22982077395351908 
H  4.569317637879383  3.4816841702181875  -0.3057213839391635 
C  2.5372178514978385  3.550909950643956  -0.10760952282878232 
H  2.4769663163441895  4.499689510292723  -0.1090058409595112 
C  1.38623306518441  2.8054969189474113  0.015988631898662505 
H  0.55394992621879  3.251435989430907  0.11455623311663164 
N  1.4109171476383717  1.445480128711564  2.670689485642905e-16 
C  2.6090807019930238  0.8942529902864789  -0.13590396933108714 
C  1.676857266007879  -2.4405052324340386  1.4888869622653509 
C  0.8003815893851329  -2.3119998062366056  2.5726347759778174 
H  0.014203366519743632  -1.786718500102851  2.4874114733038706 
C  1.0780915260637314  -2.944729351827268  3.76128248375613 
H  0.48849476181745033  -2.8429328988167675  4.499211760682307 
C  2.201185380220278  -3.724513845305364  3.886556816669584 
H  2.376004904455649  -4.171751242180709  4.706631267425297 
C  3.084677343536942  -3.8626790687719885  2.821644077625054 
H  3.8681733035309858  -4.391161598435621  2.9135170144093085 
C  2.8110683822862392  -3.2259710769296124  1.630123398221046 
H  3.405811089889033  -3.3257557554838346  0.8961541756516337 
C  1.454030558239608  -2.5734537601812146  -1.4105414033397023 
C  1.101368545108473  -3.905557662374827  -1.2726187139269514 
H  0.8973752886654792  -4.257708567039095  -0.41419212902296143 
C  1.0466746758833512  -4.725269273392368  -2.39669761926779 
H  0.7958103198329365  -5.6360965333793605  -2.3056754388500607 
C  1.355380099808396  -4.219049315424332  -3.645014071692467 
H  1.3113224474393936  -4.778709862761698  -4.411498758442154 
C  1.730213438713152  -2.887504410884129  -3.772181305044986 
H  1.963860099375199  -2.5425999392852088  -4.625258607188547 
C  1.7663255604061483  -2.064200156220473  -2.6702780974447995 
H  2.0031231171460293  -1.1500921973806864  -2.766752999402681 
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


