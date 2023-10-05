%chk=IPINULIR_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.2791171936144097  1.5633167321435533  1.6022417988769757e-17 
C  0.958942514045761  -2.8847200169196987  -0.5764583797169902 
C  0.8251468659733431  -2.7977601691775225  -2.1149102799847763 
C  2.1574066247359567  -2.4504919708415325  -2.770153039278051 
C  2.2364090179285685  -1.4137963533272315  -3.6391638956524237 
C  3.4409973616009726  -1.1306786965532774  -4.262801943490618 
C  4.519115773306639  -1.8280493429386366  -3.9583141573622598 
C  4.442433038793667  -2.829576124305054  -3.036039597139845 
C  2.4525403754299093  -1.4545112920797405  0.5208405731032654 
C  2.664589021432551  1.0501352421073231  1.0180000651309014 
C  3.463987165603025  2.029407291412955  1.6873089628347335 
C  4.563248140343754  1.6932825834790126  2.4351744254204046 
C  4.9449048070733115  0.358926910627142  2.578332754095039 
C  4.23088499669943  -0.590848739099389  1.9294807019443168 
C  3.074766892787979  -0.2665819577255881  1.146526402955071 
C  0.7384643411101942  3.168897574074638  0.6786850691616407 
C  0.5834236285933062  4.273123224262002  -0.105222853320693 
C  -0.0013000019534616936  5.428010182149567  0.42216403538761027 
C  -0.4116561914753982  5.460269479938091  1.7146465209274453 
C  -0.17450334732067296  4.4101900518563  2.5567689129992486 
C  0.4033980072703848  3.228903667575203  2.0492190310429117 
C  1.9925936477906694  1.9745820802207978  -1.6365229826699823 
C  3.369347249970321  2.186434962103462  -1.7811306080683862 
C  3.841945569820398  2.4710995670743348  -3.0440016587904073 
C  3.041399371436321  2.5573233899430203  -4.126676977906303 
C  1.659782551669966  2.4322727933617716  -3.9755415365518276 
C  1.1501034356178061  2.1089350495720605  -2.7052507225642066 
N  1.279117193614409  -1.5633167321435533  -3.608224830031759e-16 
N  3.2292914436944598  -3.144226863266652  -2.4185923988215037 
H  1.2830914366387334  -0.9307969044793564  -3.82384232324037 
H  3.666072550825361  -0.1875595066381448  -4.767991029684527 
H  5.509901274788293  -1.7070948804729733  -4.317144476682374 
H  5.2189826458692075  -3.5815720351964155  -2.73782449397852 
H  2.9363172532373927  -2.3875723634596473  0.36634258096670136 
H  3.0828869889120956  3.0099771935198003  1.475905716860827 
H  5.176949776328804  2.4976952748698977  2.802274948327366 
H  5.796544680437522  0.06322643260171268  3.2019252840308763 
H  4.402755082088383  -1.6630821359118255  2.0490031956785817 
H  1.079540285984967  4.249079900713316  -1.0680810832025815 
H  -0.1657055442691755  6.089502977395259  -0.42872641554831226 
H  -1.034451823023952  6.246395509124063  2.15262048997165 
H  -0.46527168414489206  4.373792561936276  3.6017272361050403 
H  0.553003968509141  2.338526755219491  2.6156438718697776 
H  4.067968240233426  2.142895662868386  -0.9585880308202925 
H  4.876123975775468  2.6337024067762105  -3.3344959808985797 
H  3.382011732262881  2.665428751144054  -5.1533726487512705 
H  0.9250254905723063  2.5423513823397315  -4.800731245295079 
H  0.17249175096855995  1.6209054044508937  -2.632564826269407 
H  1.6615920504884694  -3.6463007289533427  -0.23774573159887433 
H  -0.0019229090986581099  -3.212229125901066  -0.15788170557026618 
H  0.6152850140539062  -3.832366021430038  -2.3535217080300095 
H  0.04673405902519612  -2.1123625411928137  -2.372479015533233 
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


