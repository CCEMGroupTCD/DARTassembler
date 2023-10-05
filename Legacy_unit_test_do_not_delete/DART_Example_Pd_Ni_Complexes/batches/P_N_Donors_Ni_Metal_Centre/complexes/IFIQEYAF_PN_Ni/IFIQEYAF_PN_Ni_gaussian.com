%chk=IFIQEYAF_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4483704835089681  1.4079499076671722  9.380133434037627e-16 
N  1.4483704835089677  -1.4079499076671718  1.9984014443252818e-15 
N  2.8156978781159303  0.4567420431636664  0.20581526781087908 
H  3.600304567885854  0.8286028118947686  0.3535386953605007 
C  1.7768549027005065  2.465167775926635  -1.437915065426806 
C  0.7354919595177437  2.9176792832824825  -2.2462790281789102 
H  -0.1447689969042929  2.583033893958463  -2.1191815380260586 
C  0.9859399801567852  3.8626050421856144  -3.2436656578864493 
H  0.2771969958917391  4.180479589825982  -3.78944067950024 
C  2.2792307214655905  4.334185335194841  -3.431137057545925 
H  2.447801258127207  4.983717118974879  -4.102542066485504 
C  3.3276159460943493  3.8692386162560695  -2.650057945465579 
H  4.209847201975119  4.18853355533871  -2.796066585700606 
C  3.079117636263015  2.932044300737906  -1.6504979747927764 
H  3.7931124478324576  2.6101585308450677  -1.1132522207914908 
C  1.6005057354650118  2.6156912773991774  1.3454659337754231 
C  1.0261482433282236  3.8869008532892115  1.2342206849248536 
H  0.5566657980040053  4.1296152873319745  0.44359098839195066 
C  1.1430903484620512  4.792234241258949  2.2790051510584757 
H  0.7580016315805754  5.656647128658137  2.199423681929157 
C  1.8229281087673215  4.4432971645031385  3.4426284690305264 
H  1.8944133697352676  5.064083294918948  4.158301535274963 
C  2.3940284619717556  3.184400959923824  3.552465906480495 
H  2.865782280905617  2.945801260338908  4.342795899381749 
C  2.280151959133626  2.2709349814850315  2.511684321162857 
H  2.667229162348053  1.4074258513808104  2.595621148469828 
C  2.6952753957804525  -0.9066194724320742  0.1405339263675036 
C  3.825253856118953  -1.7341274637713642  0.2269446516317028 
H  4.691082730437827  -1.3606060926372432  0.3432332396436993 
C  3.654627522898105  -3.103736080170976  0.13991895257937376 
H  4.398336881708913  -3.68853018791322  0.22490072219597243 
C  2.3695972086623955  -3.616062325637455  -0.07606241166657243 
H  2.2340988129182673  -4.551703798268176  -0.17103519431604938 
C  1.3059448888503895  -2.7485409434145778  -0.14920695312165352 
H  0.43858585602346256  -3.101035646739861  -0.30976062732804366 
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


