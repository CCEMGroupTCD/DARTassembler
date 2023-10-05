%chk=EXONIFEF_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.377500640653208  1.4773598021470606  -1.2140604779440774e-15 
N  1.3775006406532082  -1.4773598021470606  -1.942890293094024e-16 
N  1.8762922373075843  -4.313326595010254  -0.4615477796330937 
C  1.3556335629959195  -2.1078002091802053  -1.354312365776319 
H  1.6427080757703445  -1.4369877005687446  -2.02393152052455 
H  0.4256589338339537  -2.3709827591825343  -1.5702601013267357 
C  2.2509532262508314  -3.3287535256157414  -1.4690516775844458 
H  2.1612055083860953  -3.7241270014278687  -2.3722865704901217 
H  3.1969905877938105  -3.0626155378157582  -1.3415999533332996 
C  2.0531460182223595  -3.74410967928136  0.8659409785205419 
H  2.9945310802764014  -3.460625774824678  0.9871265009871858 
H  1.8420053891639558  -4.424221557606302  1.5539790644737488 
C  1.130921317399822  -2.5438067859557942  1.0245106334583256 
H  0.19128168082188934  -2.8492187575988988  0.9515404960816631 
H  1.254446436094291  -2.1587559512347787  1.9285580988782454 
C  2.709731176350428  -0.8225930089838941  0.2733513605723622 
H  3.438385485697015  -1.4307310662765662  -0.008912152377793144 
H  2.805321173158152  -0.6617823251851814  1.2451530412081395 
C  2.8379150421616215  0.49468324773551164  -0.47243783218187163 
H  2.8537858495590145  0.3427985919638479  -1.4510648442828586 
H  3.6702171452197443  0.9623267723501256  -0.20957573310505398 
C  2.612802787095428  -5.560492391942228  -0.6192161182881972 
H  2.46620626134545  -5.915229358684087  -1.5217113000287106 
H  2.2982835390987617  -6.211745487905528  0.04317490492832779 
H  3.5692536403751824  -5.394765830515434  -0.48485030381855443 
C  1.1604315514566368  2.768669334345976  -1.247095502584959 
C  0.9914301429427899  4.11021185304699  -0.9127127362262936 
H  1.0630568999002068  4.389000902540162  -0.0065860764427521454 
C  0.7139942255978325  5.038262665724728  -1.9072582660313246 
H  0.6071996786825113  5.955297296036904  -1.68216872483192 
C  0.5937006359586002  4.632971230199129  -3.2261353052511965 
H  0.39746148339265597  5.27237806164747  -3.90057880836868 
C  0.7581098672133455  3.3029048393684937  -3.567993390808634 
H  0.6689775628343311  3.028120948260841  -4.4735592335301035 
C  1.0511762906796216  2.3662687518094656  -2.578002016644784 
H  1.1768107443528808  1.4538221642950884  -2.8103184366059737 
C  1.7479095907040554  2.2483412268106098  1.5904093914920325 
C  2.849949519931041  3.0886124614792996  1.748344778843228 
H  3.433909521779946  3.2479386693928465  1.0155439766380596 
C  3.095717289661615  3.697994910962234  2.975065873315752 
H  3.8382627538920606  4.282397023937864  3.0782043428623727 
C  2.255973616892959  3.4479694165206096  4.045573638899117 
H  2.4193539803355866  3.8699834574427383  4.881050404647834 
C  1.1837105775544388  2.588293071493335  3.908574309451483 
H  0.6256936403381301  2.4023709964505358  4.65572597396915 
C  0.9196347254745556  2.0003142638003313  2.685744003950248 
H  0.1707384506750349  1.4241079674310577  2.5886681901971738 
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


