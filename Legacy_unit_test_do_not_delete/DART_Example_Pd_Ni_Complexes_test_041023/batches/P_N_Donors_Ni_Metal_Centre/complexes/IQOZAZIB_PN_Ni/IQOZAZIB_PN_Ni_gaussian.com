%chk=IQOZAZIB_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  2.4621554645316  0.9722539426212742  1.050501264717623 
P  1.3644029701668061  -1.4894645128367445  -1.258990933197067e-15 
N  1.3644029701668055  1.4894645128367445  3.963781011612221e-17 
C  1.2703763242356056  2.880373381188598  -0.3180819073683682 
C  1.517753557628021  3.2982421250464924  -1.6302971897461396 
C  1.509070991588376  4.658002809201838  -1.9087111169702533 
H  1.6827511994848343  4.951262394264457  -2.794997583725738 
C  1.2536267528483083  5.590314011000982  -0.9279222182288571 
H  1.2640890021812132  6.517814433531302  -1.1377147765930002 
C  0.9866389189889988  5.182650216552086  0.3478048008730676 
H  0.7949665464779607  5.82305395009478  1.021926992823891 
C  0.9937253217687918  3.831864909040351  0.6586476484762702 
H  0.807544073685644  3.5537093192767895  1.5459604291491886 
C  1.7863573390444893  2.3380244712024885  -2.736460114389281 
C  0.7364161320394518  1.7852946458828869  -3.4576318013230996 
H  -0.15928851807763889  1.9886727407547604  -3.21111036151525 
C  0.9754863446476698  0.9394421101950554  -4.537088021516401 
H  0.24406228903288452  0.5751344577981872  -5.023430879254753 
C  2.2645140448749244  0.6261412453686451  -4.904390970246603 
H  2.4264725856743956  0.05377405054616477  -5.645088665229613 
C  3.313267558905653  1.1516426005565488  -4.189122655204017 
H  4.205863809430132  0.923233516575044  -4.426446408554614 
C  3.0831877831333223  2.0107188719092783  -3.120190923623704 
H  3.820950491058893  2.3783728074549195  -2.6462735991121806 
C  2.9312338562941362  -0.6249766993287413  0.4068209490348916 
H  3.4441579634429136  -1.1369463289694457  1.0833603038629602 
H  3.4924469786843524  -0.5231935425677472  -0.40358852786971927 
C  3.8822554179713755  2.077218223955845  1.1257329017801443 
C  4.839636770430714  2.044535481547107  0.12550555046854514 
H  4.792796605987489  1.3791007253256706  -0.5497336102793025 
C  5.867937626428552  2.972486210207061  0.10423221660637391 
H  6.525866296618393  2.9421260810396497  -0.5798255495940922 
C  5.932486072449325  3.944123105565995  1.0812077801153033 
H  6.630719860241613  4.587859765108076  1.0618795430818577 
C  4.988045150339493  3.981698542588934  2.083821928817034 
H  5.0412515313680535  4.650193203870693  2.7561648389868902 
C  3.963994397871536  3.0489358026737743  2.112764758599097 
H  3.318920904817774  3.075205555771469  2.809191164353225 
C  1.8571455093142144  0.7156034875715273  2.7173631030265533 
C  2.723489135301973  0.21161211783401912  3.692308240426879 
H  3.628964160554311  0.026649564093609728  3.4758084197384465 
C  2.2563843384567948  -0.015052940037691144  4.972064959613888 
H  2.8428246077805763  -0.34277709915198273  5.643699165151617 
C  0.9185039788541421  0.23974266719514545  5.27510304620422 
H  0.5981556444286107  0.08881232406676162  6.157361169518463 
C  0.054342849607315946  0.7144801878986714  4.3039905351453 
H  -0.8591928952889731  0.8704802322081978  4.514773027236608 
C  0.5220616842965475  0.9589295264777827  3.0269156655525924 
H  -0.0668639595982774  1.2935612217645238  2.362247007894805 
C  1.0174354684887208  -2.495794947084754  1.4660382566241341 
C  2.022509120482269  -3.282854359978931  2.0006775879897165 
H  2.8851346379311607  -3.2849053890391557  1.604380704421565 
C  1.767102525632169  -4.068019968102543  3.113728408539686 
H  2.452433410916912  -4.6137196628552495  3.476078130247387 
C  0.5197279464541303  -4.052808951119768  3.694226362181896 
H  0.33911405179714804  -4.610007318618548  4.443256415192815 
C  -0.4608478609710014  -3.236583724185977  3.1916709972623125 
H  -1.3084605212722977  -3.2033018839187055  3.6174192059554273 
C  -0.22117329817769504  -2.4584363631638624  2.064985757950016 
H  -0.9081222890659368  -1.9056505109975108  1.7115930506648358 
C  1.8156898416558285  -2.6344275786247513  -1.3182106347632856 
C  2.3402269390200865  -2.091097195104092  -2.490444787069032 
H  2.408930826316716  -1.1470013587274057  -2.5896532423867997 
C  2.759573201123752  -2.930710300500093  -3.5082943072374193 
H  3.134037676633497  -2.563551875249352  -4.300521381597187 
C  2.634980278362596  -4.299497151364436  -3.3773606580149598 
H  2.929082629811732  -4.872648668330838  -4.0764127272281385 
C  2.0873525694893496  -4.8274411461510995  -2.241794239507016 
H  1.9893876730983107  -5.76862797062273  -2.162190655943826 
C  1.6738679177699203  -4.007674713666895  -1.2064644133334357 
H  1.2917204668012854  -4.385402195692286  -0.4212919067841185 
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


