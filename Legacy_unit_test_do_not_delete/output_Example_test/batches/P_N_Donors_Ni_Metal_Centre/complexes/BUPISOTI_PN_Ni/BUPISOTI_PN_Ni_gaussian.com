%chk=BUPISOTI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.5549960940465413  -1.2892195885495994  1.5788386425153346e-16 
N  1.5549960940465413  1.2892195885495998  -1.5788386425153351e-16 
N  3.707960546514663  0.7799202916480629  0.20012139461911607 
C  2.4415229738196698  0.28842065448162013  0.150458117122782 
C  3.6002362451420122  2.141222113763899  0.08115806414339279 
H  4.322913349035394  2.7567295044736015  0.09046964242155231 
C  2.280412496061148  2.463321757368581  -0.05281563085639241 
C  1.675514742043736  3.82597525379613  -0.29980286910661486 
C  0.9050386791141067  3.7956461963203676  -1.6237665222048177 
H  0.17809608905211816  3.140805653639985  -1.5639185112433325 
H  1.5126693649027994  3.544361959812246  -2.350214746217033 
H  0.5312214687433454  4.6827083462886225  -1.8044179969409078 
C  2.8081904594130993  4.851111969634321  -0.40361579278338666 
H  3.3045027521401833  4.875613001328647  0.4396734266033676 
H  2.4303885701435792  5.737832649457452  -0.5862536542425818 
H  3.413098004128007  4.597581606845061  -1.132078816890261 
C  0.7417412561811448  4.215768315011834  0.8467587674956897 
H  1.2451731944022277  4.230683314885233  1.6868732991863196 
H  0.015423669104174031  3.5610433887875805  0.9135843141327453 
H  0.3667457860234884  5.104919089000326  0.674394593500891 
C  4.988351967091356  0.08622700390789412  0.4119718198202049 
H  4.8316219698911755  -0.8797173519337703  0.46761550225821935 
H  5.3954738465816465  0.4003056322473378  1.2470667907431148 
H  5.590638962935665  0.277467826797791  -0.3369776659912747 
C  1.9480163818645426  -1.9079232035771743  -1.7249767100211344 
C  3.444222500878685  -2.0492957219741736  -2.002780850499292 
H  3.5785390649322872  -2.3760056430310357  -2.9165941947353424 
H  3.8352180466583063  -2.688138927829103  -1.3684913806879888 
H  3.8803668132502964  -1.1780108454670526  -1.8987994927669554 
C  1.2251295034504266  -3.2234131791118985  -1.9872858316398048 
H  1.4372040365576075  -3.538703897318638  -2.890366047024946 
H  0.2586325113921166  -3.0866942046520878  -1.9068711167865358 
H  1.5156279089288092  -3.8916777809313996  -1.3319957368195403 
C  1.3973360525946719  -0.8275001923452051  -2.669897400743736 
H  1.5701058401368626  -1.0865171262840367  -3.598551199769812 
H  1.839634653763791  0.026602300961454356  -2.482138597138714 
H  0.432479613251995  -0.7339454317200648  -2.5306158936318384 
C  2.00823155828471  -2.4285211458990674  1.4071326822222294 
C  3.464618683837667  -2.913671507835207  1.4286739643770676 
H  3.6006868399535237  -3.5028702073791482  2.200354441872101 
H  4.064279724555616  -2.141141312943833  1.4962727382280618 
H  3.6566202318890273  -3.4055388209229673  0.6042491455720185 
C  1.7255853693226175  -1.6094680696206074  2.676136397781814 
H  1.9418522809132004  -2.1465372129610087  3.46699391859241 
H  0.7775833501570212  -1.3602100123787662  2.699671319876345 
H  2.2770707217848876  -0.8001717292091161  2.6722083682323796 
C  1.0812892776732914  -3.6628315877612194  1.3706783062849492 
H  1.312566738622189  -4.266109696399791  2.1073953956038545 
H  1.1941317561689573  -4.130082644408365  0.5171240929457526 
H  0.14943228827807786  -3.374985096522677  1.4649499287223395 
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
