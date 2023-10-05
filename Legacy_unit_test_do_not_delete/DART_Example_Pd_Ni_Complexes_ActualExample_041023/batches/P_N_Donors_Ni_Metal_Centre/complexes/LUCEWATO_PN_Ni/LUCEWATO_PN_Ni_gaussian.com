%chk=LUCEWATO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4045738917906747  1.4516446474602516  -2.7008006441312242e-15 
N  1.404573891790674  -1.451644647460252  1.1102230246251565e-16 
H  1.1099209791432565  -2.103875995301821  -0.6983141914361891 
C  2.845555224696955  0.4723216550962237  -0.4791623679326812 
C  4.065605605038219  1.0034274366308953  -0.9054259345784507 
H  4.189430903330707  1.945300630178723  -0.9246129817402883 
C  5.087346629467554  0.1731026962691169  -1.2974116161298057 
H  5.917482094325506  0.5401917844507955  -1.5782936399815406 
C  4.90623867052023  -1.1979607650829136  -1.28252698810866 
H  5.606557072942097  -1.769501542263864  -1.5745376304638556 
C  3.7146386548278914  -1.7417068715291277  -0.846187264582765 
H  3.601741534618167  -2.684548900667994  -0.8166805138695223 
C  2.687988449535086  -0.9060279626980897  -0.45293161324741665 
C  1.9591971642205555  2.4603843914202805  1.4100677411509976 
C  3.071865557319035  2.096799266835269  2.1496100551399513 
H  3.6106328636911336  1.3663320985408793  1.8692572673376158 
C  3.399712445538781  2.800172183744644  3.3016473668807955 
H  4.1568484066855715  2.539609737789456  3.812511965474408 
C  2.6389741954142907  3.8697671412409593  3.708166286115823 
H  2.8732366816560675  4.345904583837029  4.496261282482319 
C  1.544826312366847  4.250065519336358  2.978360391619985 
H  1.0256550820578056  4.996175813542955  3.2546380825368457 
C  1.1927343957859564  3.5455468709265676  1.8318517233923073 
H  0.4258669471032961  3.8044667075805267  1.3347206598506425 
C  1.1732488561313956  2.571995125688795  -1.418861823367715 
C  1.701913804788536  3.8592808352247903  -1.4528411418682738 
H  2.148837853684462  4.2160539862201905  -0.6940849535921572 
C  1.567139699706206  4.615871683351122  -2.611424077669902 
H  1.9158447882917198  5.4992335255509435  -2.6391505440856604 
C  0.9383116904639845  4.103191130777818  -3.7127986124243835 
H  0.8671902755098424  4.623521806278539  -4.504255226818332 
C  0.4072839276554996  2.8343463461932807  -3.677319495197037 
H  -0.03969904072949482  2.484194884615341  -4.43906227982228 
C  0.5269400787062213  2.0645957485143156  -2.523112454459274 
H  0.16182474094120924  1.1877392690535915  -2.4970723976911424 
C  1.5520982822088176  -2.2148652069397277  1.2516718027580667 
C  1.275030289730931  -3.5723268477397747  1.2440428250599078 
H  1.0190084963840413  -4.011051921478895  0.4414179556445217 
C  1.3783179514849353  -4.269950124050396  2.4221262572832054 
H  1.1709324037096192  -5.196393236249332  2.435070925416657 
C  1.773234849466577  -3.6576886855812853  3.5735060836881716 
H  1.8544123288163463  -4.159686489434279  4.375647799205576 
C  2.0562299403435556  -2.2994433474748823  3.5707420193340975 
H  2.3257782149732984  -1.8680706343969424  4.372736732568482 
C  1.944398881419944  -1.5696530252679417  2.3898405162978973 
H  2.1377956383350982  -0.6395282923573892  2.3758521052811505 
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


