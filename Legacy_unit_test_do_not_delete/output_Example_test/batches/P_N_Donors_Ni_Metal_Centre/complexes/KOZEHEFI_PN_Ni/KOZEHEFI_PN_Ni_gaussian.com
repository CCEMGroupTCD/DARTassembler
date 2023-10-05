%chk=KOZEHEFI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3122583453725891  -1.5356034758361283  -8.373191619686828e-16 
N  1.3122583453725873  1.5356034758361283  3.670543241683518e-16 
N  2.61033773106704  -1.0433422300816495  -1.0020896588446664 
C  2.4986014359807367  1.4113331951270451  -0.4963632829882142 
H  3.031825572221438  2.198640686301946  -0.4864358024294062 
C  3.1162776886748826  0.25903086067411907  -1.0561131886354302 
C  4.294654651995481  0.21733590100027644  -1.7680535742388472 
H  4.855309928291511  0.9589777590195754  -1.958126574938878 
C  4.522247608057347  -1.1070075489861713  -2.1645719380368154 
H  5.255039351949386  -1.4206351798026267  -2.683042003973256 
C  3.507894313359798  -1.857364785168256  -1.6746860618927657 
H  3.424046110340498  -2.7973434565593567  -1.7769652688717827 
C  1.001379106642433  2.838237379364401  0.549074159203711 
C  0.8111362333645965  3.924058567084827  -0.3217916498991665 
C  0.5425764138160863  5.158885245871501  0.2560753388454322 
H  0.42398143761880325  5.915033947845336  -0.3056768326003549 
C  0.4440146917281155  5.3122952919589075  1.6156459308701934 
H  0.26730536891208123  6.169562412485066  1.986488268386663 
C  0.6039496582653858  4.216365371584337  2.447200490779845 
H  0.5284442906580877  4.329025181048367  3.3863614816616128 
C  0.8711122863731738  2.9589955489817434  1.9350019459255359 
C  0.8309754308841787  3.7886172710257946  -1.8314074955909359 
H  0.95079941131822  2.819824816847502  -2.0544411128779765 
C  2.008063294461352  4.561458739331151  -2.435311657099188 
H  1.9130841197195667  5.515940018557176  -2.2288154489683336 
H  2.015646420686295  4.436504818596163  -3.407159250032174 
H  2.847305357105572  4.226217267101677  -2.056510765355103 
C  -0.489451623740091  4.253087963638445  -2.4453195569810413 
H  -1.2288810953850597  3.750138717851234  -2.0468428040519364 
H  -0.47149737945339854  4.097469567509274  -3.412625139461728 
H  -0.6129527431475688  5.208017193791143  -2.2706174994514314 
C  1.0323119301860146  1.756179663198279  2.8504030382023333 
H  0.7669754618896888  0.9475068452861904  2.3279599674643316 
C  2.4798251005295358  1.5758244677600817  3.253885720805391 
H  3.0331406721018372  1.4749685127970529  2.451136342353965 
H  2.5663266336660655  0.7761243209944277  3.8122611232256802 
H  2.77855959516447  2.361502264134613  3.759848171980659 
C  0.1323817016478046  1.8171372290497305  4.079817739972545 
H  0.3940574439599651  2.579181124764953  4.639285028460939 
H  0.22470188346702735  0.988457732666644  4.594515910071423 
H  -0.8000589138916618  1.9238808758174935  3.796768994405439 
C  2.0556979598582004  -1.9123113917159251  1.6040472485955692 
C  3.4253125236717263  -2.0940947618323587  1.7501222051728182 
H  3.999470949158528  -2.0327151632034353  0.9957114302532872 
C  3.946125902680408  -2.366000131004991  3.0049954634984175 
H  4.881809093961905  -2.4920890224713075  3.108863442212293 
C  3.1216789311566786  -2.4569009020676478  4.101120730502885 
H  3.489661169287719  -2.6451679837389657  4.956810410456436 
C  1.7621841125670001  -2.274860597132945  3.9600602271463963 
H  1.1947633315459858  -2.3368307225575338  4.720341728184266 
C  1.2231885743992443  -2.000441988012155  2.715763013132188 
H  0.2864785885191643  -1.8738656041613921  2.6197434765584946 
C  0.8373449146335709  -3.0863006821597962  -0.7772071651777821 
C  0.9492427522699427  -4.286740143913606  -0.10241832220796784 
H  1.2810323028437347  -4.3057908789896455  0.7858976877924846 
C  0.5724374848030891  -5.469295757113969  -0.7323388426687297 
H  0.6449264882410065  -6.2972026810076915  -0.2732961978639008 
C  0.09577681513316261  -5.434889107981297  -2.0204983780079435 
H  -0.16303773832142143  -6.243492217142601  -2.4486994741938077 
C  -0.011753217186629916  -4.231721139118765  -2.702592529025612 
H  -0.3324521203449582  -4.218519749359627  -3.5962893102783497 
C  0.3460195594253024  -3.0615202941242217  -2.0823184584008314 
H  0.2589543921366424  -2.2349794422881657  -2.5417064936075935 
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

-Ni 0
lanl2dz