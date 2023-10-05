%chk=ULEKIHUJ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
C  2.1002984326040264  -0.4549634097392712  1.4574064047854909 
H  1.477782273813519  -0.5106997144112844  2.198171304322612 
H  2.954801690574434  -0.8146805443769842  1.7423929398463425 
C  2.2569169179789124  0.9812608555998867  1.0047064603522293 
C  3.192950740994215  1.8446576274904096  1.6251740956567593 
H  3.681071551358375  1.539273098307401  2.3553689078987294 
C  3.3908300296084515  3.122277272203495  1.1680536951267286 
C  2.6794227879804415  3.5497915564970537  0.004079892464394885 
C  2.6341737044227767  -1.3558370041932553  -2.512560130458381 
H  1.8287982267477405  -1.0012871445495821  -2.814803851627829 
C  1.7126275039998191  2.647716636464134  -0.5771699742198528 
C  1.0458038293774732  3.027537728871856  -1.7489338346991596 
H  0.4121062885097355  2.453974230378006  -2.116895517342444 
C  1.304954871511728  4.230912769624108  -2.3714220102135575 
H  0.8850258979714365  4.448519165228918  -3.1726824446375153 
C  2.2146090134333756  5.126042625974582  -1.7738299902243508 
H  2.3518011590574064  5.95934801569899  -2.165212019238603 
C  2.899131593714699  4.810352391551486  -0.6439688956755922 
H  3.5136589016954387  5.412346793350032  -0.2899156963549795 
C  1.0675210594955207  -3.0610918191760437  0.6332049602819638 
C  1.4488382127082904  -4.207163423881093  -0.04981412843702553 
H  1.9118647973092622  -4.1333909144371175  -0.8523177600818443 
C  1.1424018099633966  -5.471014487472494  0.46174990482925554 
H  1.4021204919977557  -6.233104342952263  -0.0039047938232571675 
C  0.4572151270121516  -5.595889290393136  1.6497014918486312 
H  0.2634993280402522  -6.438305430463573  1.9955080796654447 
C  0.060275795586463454  -4.458009647652647  2.3216062066736494 
H  -0.4072523818222826  -4.53604371342797  3.121061100903114 
C  0.3489130090359225  -3.2067721903516047  1.82006543618838 
H  0.06158112397929716  -2.4507854081544105  2.2795969399755918 
C  2.8433162223124935  -1.626359304127899  -1.1626598258251761 
C  4.06269920732956  -2.1403810648104153  -0.7385344134815155 
H  4.2030508435893745  -2.3186308932222737  0.16317517772307466 
C  5.062083232597711  -2.3845790274344996  -1.6480507767698112 
H  5.883508418998477  -2.703229178414043  -1.3511200366515452 
C  4.8673113505396275  -2.166792764375532  -2.994104623234502 
H  5.534118094123003  -2.376121370436014  -3.607967524374911 
C  3.6895172970489174  -1.6408757037216546  -3.4094293189547815 
H  3.57213251190547  -1.4615861938141623  -4.314149360654961 
C  4.36710948778318  4.048887894633731  1.8693158101347447 
H  4.367636796220802  4.902417571338556  1.4308331389492468 
H  5.249082832534691  3.6706554333242734  1.8352276250899844 
H  4.100858689442744  4.159248746302258  2.78492883691085 
N  1.4551439533599415  1.4009482770609338  4.8277516686117754e-17 
P  1.4551439533599415  -1.4009482770609338  -5.551115123125783e-17 
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

