%chk=CUBISASO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
C  2.8396012509190696  0.4319071734416376  0.4808581174897794 
H  3.669298911753593  0.8219807256699415  0.10786789582473845 
H  2.9224309648048568  0.4086186815404202  1.4674147983969594 
C  2.6473289154404647  -0.9591401262521897  -0.04756661605712797 
C  3.713041677398383  -1.6797348986287326  -0.5568395026549575 
H  4.582735359343072  -1.298652772519552  -0.5828487845672288 
C  3.496946134169993  -2.9596498451624815  -1.0256374397625496 
H  4.215827543336215  -3.4727295005586787  -1.3731855325119129 
C  2.2214261249111473  -3.4825471946376916  -0.9832716765550624 
H  2.0484841592691954  -4.359837739720291  -1.3038215104788553 
C  1.2079447541735184  -2.711763116284106  -0.4687115093011024 
H  0.331628481849807  -3.0787400216737937  -0.4421242610176969 
C  1.7437854157940478  1.8535553390081485  -1.7793262685351219 
H  1.6826659721230222  1.0126788562477076  -2.298042698143612 
H  1.030364097535779  2.4560702739521107  -2.1088148204214123 
C  3.0670643900015793  2.4909420275579293  -2.072016865718994 
N  4.095208862306845  1.6592390754636728  -2.2866947489437206 
C  5.278035484098654  2.215307695894803  -2.5748160465813803 
H  6.018372667053611  1.6386430494810744  -2.7257901871182737 
C  5.487442553098713  3.577055795028477  -2.6659809738210947 
H  6.342964420599643  3.9233023524526356  -2.889304048535013 
C  4.426260822229803  4.423428715956293  -2.426241734666936 
H  4.538222358845428  5.366111483642802  -2.4710248098803977 
C  3.1977141156246653  3.874585167389982  -2.1191688623925935 
H  2.4517968393533778  4.434343729728948  -1.9433500669006034 
C  1.6910435484735806  3.025986282151485  0.8783283326107318 
C  0.5919657658941462  3.8202834314837264  1.1776389010893156 
H  -0.27390355524097454  3.5520576750635415  0.8956080653464088 
C  0.7502261517857594  5.004625899913725  1.8860396053030941 
H  -0.008127016688243671  5.542227139020646  2.084638471302312 
C  2.0039824707393  5.402624637241924  2.3033208027817906 
H  2.1083914340586447  6.206164489265223  2.7985067522886045 
C  3.1046608099915636  4.6291717773342285  1.9994710915732659 
H  3.968513395622508  4.908586013868957  2.2782629937520755 
C  2.9566330112571197  3.445923228993782  1.2892990887091202 
H  3.71934631195179  2.919474948584825  1.0825343939522347 
N  1.3938601911597879  -1.4619349395578443  2.9185534355668875e-15 
P  1.3938601911597877  1.4619349395578447  -1.7903539442911937e-16 
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


