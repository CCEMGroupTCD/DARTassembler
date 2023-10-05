%chk=OHODAVEQ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4045738917906718  -1.4516446474602525  -1.2023034095144315e-15 
N  1.404573891790674  1.451644647460252  -2.887974995636742e-16 
H  1.4326479923189153  2.103875995301821  0.7574132912567028 
C  2.9130492618053507  -0.47232165509622437  -0.17471644080133947 
C  4.198937189664895  -1.0034274366308962  -0.3040060219961689 
H  4.319269818258982  -1.9453006301787243  -0.3389474840670549 
C  5.290609343804491  -0.17310269626911762  -0.38055276209279654 
H  6.16167345257495  -0.5401917844507963  -0.4768176038049328 
C  5.120179274695255  1.1979607650829136  -0.31750328557777624 
H  5.878292326297208  1.769501542263864  -0.34881911241024877 
C  3.8558177659553694  1.7417068715291277  -0.20936944748063124 
H  3.741028135079228  2.684548900667995  -0.1883992607236249 
C  2.759159670548663  0.9060279626980895  -0.1318989814054473 
C  1.3113129047637573  -2.460384391420282  -1.5123492973695825 
C  2.0071888465967525  -2.0967992668352706  -2.6528362375659795 
H  2.6139600596177344  -1.3663320985408804  -2.6264432253649788 
C  1.8174470281976998  -2.800172183744646  -3.8354907061699977 
H  2.287744536988021  -2.539609737789458  -4.618470753797963 
C  0.9561817092014291  -3.8697671412409615  -3.8824200913431466 
H  0.8354322872419544  -4.345904583837032  -4.695680328201281 
C  0.2729762611658073  -4.25006551933636  -2.7585844497155505 
H  -0.3143126651616923  -4.996175813542957  -2.789565829479011 
C  0.4384081158529112  -3.5455468709265694  -1.5706942420374115 
H  -0.04651315858182614  -3.8044667075805285  -0.7960482998239996 
C  1.794559128078045  -2.571995125688797  1.383687703714607 
C  2.288052568299177  -3.859280835224792  1.1910599630251477 
H  2.3724390926532544  -4.216053986220192  0.3145150564537468 
C  2.6555440537830917  -4.6158716833511235  2.2980506976874966 
H  2.983295901988287  -5.499233525550945  2.175810271708963 
C  2.5510933237035225  -4.1031911307778195  3.5619892151791666 
H  2.8211194501238746  -4.62352180627854  4.309349716666534 
C  2.05482460428683  -2.834346346193282  3.754256245013017 
H  1.9716468708985853  -2.484194884615342  4.633532827564416 
C  1.6754809325054574  -2.0645957485143165  2.657620541576208 
H  1.333569055195535  -1.187739269053592  2.7883246447335006 
C  1.009297034053013  2.214865206939728  -1.1967464030967092 
C  0.7614123002253446  3.5723268477397756  -1.0727382078367498 
H  0.8685816824061557  4.011051921478895  -0.23711355334431145 
C  0.3571431400660796  4.269950124050396  -2.184095648309325 
H  0.1637173500530753  5.196393236249333  -2.1081825821283307 
C  0.22846525910062132  3.657688685581286  -3.3944992437740837 
H  -0.03696269671433483  4.159686489434278  -4.155793612069749 
C  0.48611405772822325  2.2994433474748823  -3.5115930441137 
H  0.39147014637510735  1.868070634396942  -4.352363121159487 
C  0.8838312386933479  1.5696530252679415  -2.3940709684577897 
H  1.065019983421595  0.6395282923573885  -2.4631261638588056 
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