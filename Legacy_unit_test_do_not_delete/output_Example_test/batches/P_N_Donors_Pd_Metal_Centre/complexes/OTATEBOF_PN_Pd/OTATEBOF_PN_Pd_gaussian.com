%chk=OTATEBOF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.682264204130849  1.2892195885495994  2.245560257031743e-16 
N  1.682264204130849  -1.2892195885495998  -2.636779683484747e-16 
N  3.835228656598971  -0.7799202916480629  -0.20012139461911532 
C  2.5687910839039754  -0.28842065448162013  -0.15045811712278384 
C  3.727504355226319  -2.141222113763899  -0.0811580641433937 
H  4.450181459119701  -2.7567295044736015  -0.09046964242155635 
C  2.4076806061454548  -2.463321757368581  0.052815630856391446 
C  1.8027828521280436  -3.82597525379613  0.2998028691066147 
C  1.0323067891984148  -3.795646196320368  1.6237665222048154 
H  0.305364199136426  -3.140805653639985  1.5639185112433314 
H  1.6399374749871087  -3.544361959812247  2.3502147462170337 
H  0.6584895788276524  -4.6827083462886225  1.8044179969409067 
C  2.935458569497408  -4.851111969634321  0.4036157927833845 
H  3.431770862224491  -4.875613001328647  -0.4396734266033685 
H  2.5576566802278853  -5.737832649457452  0.5862536542425767 
H  3.540366114212315  -4.597581606845061  1.1320788168902622 
C  0.8690093662654522  -4.215768315011834  -0.846758767495692 
H  1.3724413044865365  -4.230683314885233  -1.6868732991863211 
H  0.1426917791884823  -3.5610433887875805  -0.9135843141327464 
H  0.49401389610779467  -5.104919089000326  -0.6743945935008936 
C  5.115620077175663  -0.08622700390789384  -0.4119718198202033 
H  4.9588900799754825  0.8797173519337709  -0.4676155022582202 
H  5.522741956665955  -0.4003056322473372  -1.247066790743115 
H  5.717907073019973  -0.27746782679779103  0.33697766599127543 
C  2.0752844919488496  1.907923203577174  1.724976710021134 
C  3.571490610962992  2.0492957219741728  2.0027808504992923 
H  3.7058071750165933  2.376005643031035  2.916594194735338 
H  3.9624861567426133  2.688138927829103  1.368491380687988 
H  4.007634923334605  1.178010845467052  1.8987994927669567 
C  1.3523976135347338  3.223413179111897  1.9872858316398059 
H  1.56447214664191  3.5387038973186367  2.8903660470249455 
H  0.38590062147642556  3.0866942046520873  1.9068711167865373 
H  1.6428960190131177  3.8916777809313996  1.3319957368195416 
C  1.5246041626789797  0.8275001923452043  2.6698974007437353 
H  1.6973739502211718  1.0865171262840359  3.598551199769814 
H  1.9669027638480998  -0.026602300961455105  2.4821385971387113 
H  0.5597477233363013  0.7339454317200635  2.530615893631839 
C  2.13549966836902  2.428521145899068  -1.407132682222229 
C  3.5918867939219714  2.913671507835208  -1.4286739643770656 
H  3.7279549500378293  3.5028702073791487  -2.200354441872099 
H  4.191547834639927  2.1411413129438333  -1.4962727382280574 
H  3.7838883419733333  3.4055388209229673  -0.604249145572018 
C  1.8528534794069245  1.6094680696206085  -2.6761363977818142 
H  2.069120390997509  2.1465372129610105  -3.466993918592411 
H  0.9048514602413296  1.3602100123787668  -2.699671319876348 
H  2.404338831869192  0.8001717292091168  -2.6722083682323787 
C  1.2085573877575966  3.6628315877612203  -1.3706783062849508 
H  1.4398348487064958  4.266109696399791  -2.107395395603854 
H  1.3213998662532658  4.130082644408365  -0.517124092945751 
H  0.27670039836238614  3.374985096522678  -1.4649499287223375 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
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

-Pd 0
lanl2dz
