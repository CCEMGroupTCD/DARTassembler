%chk=ZIKUXUHI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  2.4621554645316  -0.9722539426212741  -1.0505012647176233 
P  1.3644029701668061  1.4894645128367445  1.4413977280059761e-15 
N  1.3644029701668055  -1.4894645128367445  -2.220446049250313e-16 
C  1.2703763242356056  -2.880373381188598  0.3180819073683679 
C  1.517753557628021  -3.2982421250464924  1.6302971897461391 
C  1.509070991588376  -4.658002809201838  1.9087111169702526 
H  1.6827511994848343  -4.951262394264457  2.7949975837257375 
C  1.2536267528483083  -5.590314011000982  0.9279222182288565 
H  1.2640890021812132  -6.517814433531302  1.1377147765929994 
C  0.9866389189889988  -5.182650216552086  -0.3478048008730682 
H  0.7949665464779607  -5.82305395009478  -1.0219269928238917 
C  0.9937253217687918  -3.831864909040351  -0.6586476484762707 
H  0.807544073685644  -3.5537093192767895  -1.545960429149189 
C  1.7863573390444893  -2.338024471202489  2.7364601143892804 
C  0.7364161320394518  -1.7852946458828873  3.4576318013230996 
H  -0.15928851807763889  -1.9886727407547609  3.2111103615152494 
C  0.9754863446476698  -0.939442110195056  4.537088021516401 
H  0.24406228903288452  -0.5751344577981878  5.023430879254753 
C  2.2645140448749244  -0.6261412453686457  4.904390970246603 
H  2.4264725856743956  -0.053774050546165464  5.645088665229613 
C  3.313267558905653  -1.1516426005565492  4.189122655204017 
H  4.205863809430132  -0.9232335165750446  4.426446408554614 
C  3.0831877831333223  -2.010718871909279  3.1201909236237038 
H  3.820950491058893  -2.37837280745492  2.64627359911218 
C  2.9312338562941362  0.6249766993287413  -0.4068209490348915 
H  3.4441579634429136  1.136946328969446  -1.08336030386296 
H  3.4924469786843524  0.5231935425677472  0.4035885278697193 
C  3.8822554179713755  -2.077218223955845  -1.1257329017801445 
C  4.839636770430714  -2.044535481547107  -0.1255055504685454 
H  4.792796605987489  -1.3791007253256706  0.5497336102793022 
C  5.867937626428552  -2.972486210207061  -0.10423221660637427 
H  6.525866296618393  -2.9421260810396497  0.5798255495940918 
C  5.932486072449325  -3.944123105565995  -1.0812077801153037 
H  6.630719860241613  -4.587859765108076  -1.0618795430818584 
C  4.988045150339493  -3.9816985425889335  -2.0838219288170343 
H  5.0412515313680535  -4.650193203870693  -2.7561648389868907 
C  3.963994397871536  -3.048935802673774  -2.1127647585990976 
H  3.318920904817774  -3.0752055557714684  -2.8091911643532255 
C  1.8571455093142144  -0.7156034875715269  -2.7173631030265533 
C  2.723489135301973  -0.21161211783401868  -3.692308240426879 
H  3.628964160554311  -0.0266495640936093  -3.4758084197384465 
C  2.2563843384567948  0.015052940037691753  -4.972064959613888 
H  2.8428246077805763  0.3427770991519834  -5.643699165151617 
C  0.9185039788541421  -0.2397426671951448  -5.27510304620422 
H  0.5981556444286107  -0.08881232406676087  -6.157361169518463 
C  0.054342849607315946  -0.7144801878986708  -4.3039905351453 
H  -0.8591928952889731  -0.8704802322081973  -4.514773027236608 
C  0.5220616842965475  -0.9589295264777824  -3.0269156655525924 
H  -0.0668639595982774  -1.2935612217645236  -2.362247007894805 
C  1.0174354684887208  2.495794947084754  -1.466038256624134 
C  2.022509120482269  3.2828543599789315  -2.000677587989716 
H  2.8851346379311607  3.2849053890391557  -1.6043807044215646 
C  1.767102525632169  4.068019968102543  -3.1137284085396857 
H  2.452433410916912  4.6137196628552495  -3.4760781302473864 
C  0.5197279464541303  4.0528089511197685  -3.6942263621818956 
H  0.33911405179714804  4.610007318618549  -4.443256415192814 
C  -0.4608478609710014  3.2365837241859774  -3.191670997262312 
H  -1.3084605212722977  3.203301883918706  -3.617419205955427 
C  -0.22117329817769504  2.458436363163863  -2.0649857579500157 
H  -0.9081222890659368  1.905650510997511  -1.7115930506648356 
C  1.8156898416558285  2.6344275786247513  1.3182106347632858 
C  2.3402269390200865  2.0910971951040915  2.4904447870690323 
H  2.408930826316716  1.1470013587274055  2.5896532423867997 
C  2.759573201123752  2.9307103005000927  3.5082943072374198 
H  3.134037676633497  2.5635518752493516  4.300521381597187 
C  2.634980278362596  4.299497151364436  3.37736065801496 
H  2.929082629811732  4.872648668330837  4.076412727228139 
C  2.0873525694893496  4.8274411461510995  2.2417942395070165 
H  1.9893876730983107  5.76862797062273  2.162190655943827 
C  1.6738679177699203  4.007674713666895  1.2064644133334361 
H  1.2917204668012854  4.385402195692286  0.42129190678411904 
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


