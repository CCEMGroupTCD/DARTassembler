%chk=TIGODAWO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4274310000837163  1.4291748458463713  -1.9322549502024548e-17 
N  1.4274310000837165  -1.4291748458463718  1.3877787807814457e-16 
C  2.968441449004998  0.45129836413964686  0.18167386394402876 
H  3.6704097642874656  0.8207385145414801  -0.41112188470013356 
H  3.2902537229922633  0.5046573847112961  1.116245688489792 
C  2.6964793540077627  -0.9780961749100531  -0.1817793113339092 
C  3.693366350945431  -1.8161877656493668  -0.6685015225269842 
H  4.574566807891624  -1.4857068161891278  -0.7940576323811211 
C  3.395713028152098  -3.1323065564793806  -0.9686719538726241 
H  4.067241451810849  -3.711332609453277  -1.3104030414193473 
C  2.1126580944734763  -3.5963997511837738  -0.7639438322281776 
H  1.8873565909710215  -4.499256482639782  -0.9558275860587305 
C  1.1610503057715313  -2.7233284166116514  -0.27707731948544334 
H  0.28063045279734244  -3.0481735965763592  -0.13054363208688743 
C  1.3066089937232204  1.859933287886693  -1.7650203447524644 
C  0.3049733373451087  2.752819000022416  -2.153193729898706 
H  -0.2952102455112988  3.109010855399306  -1.5087961359919426 
C  0.19169101866314797  3.117528945425437  -3.4927965747740926 
H  -0.47899946727512566  3.736256984594307  -3.758241664421889 
C  1.0507339563191906  2.5836014471786437  -4.441552453743611 
H  0.9644415078908537  2.831230371679069  -5.355559291297376 
C  2.0273375772632276  1.6962564706799999  -4.058110606468635 
H  2.6141014314903805  1.3308784124825417  -4.709952967749633 
C  2.163698781425122  1.3259534937292785  -2.720384738713599 
H  2.8399993694310757  0.7107646856171952  -2.4636878259225647 
C  1.7750388274613134  3.030474278628303  0.7834208415137486 
C  0.7881727947776569  3.673414967565376  1.5198093015295078 
H  -0.04933541304353062  3.2488036670210945  1.6603695565173642 
C  1.0214232977027649  4.939398290220698  2.0518866357994203 
H  0.34375724530814256  5.380253082642456  2.55170520631344 
C  2.240781922265903  5.549996806235866  1.850414436864661 
H  2.40166678011476  6.413956393537985  2.212129012218913 
C  3.2264824567815573  4.917983289747557  1.1279769093014371 
H  4.065902668342032  5.343859341071111  1.0005511374464364 
C  3.0019988621081404  3.657399201498318  0.5816065841781476 
H  3.680899666432607  3.228490351205067  0.07447980907379984 
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