%chk=XEBIKARU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4516608040448018  -1.404557193566713  1.7200864713208645e-16 
N  1.4516608040448018  1.4045571935667125  -1.720086471320864e-16 
N  2.6788019702865933  1.2693458964619362  0.6096324493086152 
C  3.605274006884248  2.05108425858001  -0.010487457318788836 
C  5.046053382868854  2.0803802674491907  0.37968711079373446 
H  5.138044312580941  2.5105572675510395  1.2323641524889404 
H  5.546469703362158  2.566112627174273  -0.2796288839309957 
H  5.38034067643777  1.1820808847120419  0.4374213669371842 
C  1.6151282780537892  2.273417540151123  -1.0008603264563751 
C  2.949777734080465  2.7017545511006356  -1.028808932478526 
H  3.3219548615405587  3.311740711809658  -1.6240325539343734 
C  0.5002010932989417  2.643940187704158  -1.9070871801442113 
H  -0.21090899208860558  2.0054921101838405  -1.8176165214300448 
H  0.813515522413473  2.6485736052892936  -2.8147951838271976 
H  0.17630408759596228  3.5180222013261178  -1.676696635287816 
C  2.9196269445396106  0.29258120771603946  1.6273144647707989 
C  2.5457715816649635  -1.0492321400296205  1.4066317117939378 
C  2.890896841126044  -2.000703576095798  2.359755057238856 
H  2.6675248942277743  -2.89266711146357  2.2223356084220574 
C  3.566554960102164  -1.6283575994453567  3.5151273002557444 
H  3.8193203560162505  -2.2744795901743133  4.135081151542147 
C  3.864667455397015  -0.30203009314469065  3.7453952879990373 
H  4.284653391263949  -0.055021540136613534  4.5381452412723275 
C  3.5433292071879414  0.6746941771687912  2.801797826144399 
H  3.7452646367636646  1.5690003796882261  2.959152611318964 
C  2.40641182583689  -1.1015564528711466  -1.4991712468977567 
C  3.776227599622979  -1.4109724538646  -1.5551116993100182 
H  4.2122303127296075  -1.7396386926723921  -0.8014637525659115 
C  4.467703724584661  -1.2278474724269648  -2.7280596008732125 
H  5.376346259815556  -1.4252695220736955  -2.761612804632483 
C  3.8305377047920848  -0.7531040778153675  -3.8535563456461044 
H  4.309828573567964  -0.637574563785441  -4.641384045397806 
C  2.4670410962203486  -0.4431295805188855  -3.8178657384893704 
H  2.0378908918779723  -0.12412062022002612  -4.577838248098747 
C  1.7627920345923056  -0.6171664735123706  -2.633665337766065 
H  0.8566935135652634  -0.4097754528161647  -2.5982907038145004 
C  1.0799383778404  -3.1776873881982004  0.010041949154549396 
C  1.448214047699861  -3.975304434432629  -1.0846183353839298 
H  1.9225127995386222  -3.6081324514317776  -1.795680244716325 
C  1.0925136470762604  -5.326761739762958  -1.0940309261803147 
H  1.3490280810832351  -5.8659662545159765  -1.8078722645854628 
C  0.3672022422752268  -5.868376565762473  -0.05760633831918999 
H  0.11921757537918576  -6.764741353018251  -0.0828496176696924 
C  0.005242067443747667  -5.0785100618696895  1.024815334534872 
H  -0.47621683686702765  -5.45088282276229  1.7287685079811645 
C  0.35723912119643186  -3.738348658246081  1.0653352850164735 
H  0.11221258462511163  -3.2144890437212057  1.792961116502493 
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

