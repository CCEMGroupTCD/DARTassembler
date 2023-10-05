%chk=CANOKAYI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.300791176553715  1.5453291930847608  4.6832680838715325e-15 
N  1.300791176553716  -1.5453291930847608  2.220446049250313e-15 
C  0.6815512784702766  -2.237668178064591  1.1637563226460914 
C  1.720515405355544  -2.6367997208841305  -0.9263802290107965 
C  2.504344590991609  -0.8518650594907877  0.48368136620474234 
C  2.180710374709876  0.4171414756704459  1.1597194182179391 
C  2.721992752416646  2.2278831141004414  -0.9686707519435496 
C  0.6127947287213404  2.8770587230834317  1.099951592889251 
C  2.3364830648240917  3.020036064732532  -2.1947498821701443 
C  3.752847519986482  3.00434253476589  -0.1462887809358866 
C  0.3503283767381564  4.187524751333067  0.3688328864671823 
C  1.4267788697771788  3.1147522589615524  2.3755955993077342 
H  0.39016011398939654  -1.5926325769379348  1.7952937107251064 
H  -0.06242404855030692  -2.7491491968255217  0.8642384096059652 
H  1.320993683578152  -2.8120880259340613  1.5655470379737118 
H  2.126784311115962  -2.25532176289433  -1.6933215590877306 
H  2.340780439023396  -3.2016622006108832  -0.48196715876423585 
H  0.9573513429855214  -3.1386254862098912  -1.1833772234720439 
H  3.0691554009120416  -0.6734988471436403  -0.25391454843279304 
H  2.9534541094952775  -1.4191450500085852  1.1000087327458699 
H  1.6281106510636048  0.2425145992692319  1.9077382277838726 
H  2.9877410225281595  0.8298277743916658  1.4438867673800055 
H  3.136946863383526  1.4318039482469285  -1.270534392492825 
H  -0.22409843892369352  2.529987673497312  1.3728047630893954 
H  1.6994978377085033  2.5386517691863877  -2.6959686771470044 
H  1.9701429857619086  3.8559185716963076  -1.9251851570582166 
H  3.1124699471012676  3.176110787244276  -2.7215909013843427 
H  3.986920089929545  2.4908859369415066  0.6214850882149665 
H  4.528829304060551  3.1507528245685177  -0.6748643647696049 
H  3.3838665688935965  3.8311462576991944  0.12250039569558413 
H  -0.1542205868267521  4.0067741279384315  -0.4146060553196296 
H  -0.1416823538750993  4.764909300847016  0.9392263059102819 
H  1.1747866899479111  4.587196049697081  0.13738597453408746 
H  1.5766941316434038  2.291738984061075  2.8071517991209016 
H  2.2541575935138187  3.5201198064091863  2.1431075587890027 
H  0.938232019256126  3.698344408152597  2.946859214680331 
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
