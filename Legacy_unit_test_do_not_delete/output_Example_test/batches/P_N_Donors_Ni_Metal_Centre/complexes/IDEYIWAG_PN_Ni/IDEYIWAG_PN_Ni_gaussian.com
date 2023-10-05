%chk=IDEYIWAG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2008831822871033  -1.6241858214194578  -4.424449332108153e-16 
N  3.7365077909015882  -1.721162416044785  -1.2873241061534326 
N  1.2008831822871033  1.6241858214194582  -4.209500016672167e-16 
C  2.8091120337966142  -1.0109226271074296  -0.5831368498686044 
C  3.1992790307885715  0.31415204514842715  -0.4970225223480863 
C  4.441639681100051  0.4157628047201168  -1.230644671752605 
C  4.734221975775226  -0.8447097058764901  -1.7205784334724161 
C  5.850221286449177  -1.1339316987082348  -2.5153946307391974 
C  6.689340950943956  -0.06910212244182147  -2.7771815138454246 
C  6.426254278470786  1.2046515070486523  -2.2948720653304164 
C  5.313402280469411  1.4754653506368136  -1.5577881334302788 
C  2.5434306405783267  1.4231599531387886  0.1523224622891161 
C  0.643209258885595  2.7690741774665995  0.4574194326279279 
C  1.3106768332549168  3.770591933962212  1.082045897476759 
C  2.6957821635787695  3.5529098342065284  1.3466797905336583 
C  3.4438355203416866  4.542741877312885  2.0714081823299675 
C  4.7481459306702725  4.250847936516339  2.435741190489692 
C  5.313170439934561  3.021993479915007  2.1256670083627367 
C  4.623904436000691  2.085637032922456  1.3887834626131292 
C  3.289479664975613  2.3725822142273216  0.9592672889937851 
C  0.8159253557488806  -3.119832152623723  -0.9598579133091582 
C  0.8186166934398327  -4.388730161595907  -0.40021965556772193 
C  0.42506410987257803  -5.47491378463616  -1.1298699217958013 
C  0.039388239851329  -5.329382462491368  -2.449139040496824 
C  0.05067680137213548  -4.0667355793855124  -3.029450126865058 
C  0.4143073452356605  -2.954381170072263  -2.295773511692594 
C  1.4753372933889581  -2.2081978737910397  1.711059202244065 
C  2.7306974698840696  -2.1674602906109475  2.323627187521502 
C  2.860000304036303  -2.510338588559207  3.6523121237935725 
C  1.7493050319704315  -2.917203943055697  4.383691407011951 
C  0.5093602151456809  -2.9507961923417523  3.7788649260958844 
C  0.35850929051843083  -2.6191852073811206  2.4573233402225605 
C  3.8667498494030865  -3.1915768129102915  -1.4596790239478985 
H  6.015297004327818  -1.9874638264474414  -2.84561462450675 
H  7.451456361661924  -0.20945935065648855  -3.291308190097622 
H  7.024907439569619  1.8910760345239934  -2.4794377739186735 
H  5.13305334285804  2.3427945603998097  -1.2758665794576616 
H  -0.27170706290464386  2.8773287704332704  0.33405696984665867 
H  0.8862135661017593  4.560661621849604  1.3261847782786256 
H  3.0608828304405074  5.362062938025635  2.2920592117179805 
H  5.250542765367577  4.885593097078698  2.8945881820116823 
H  6.1742905270179715  2.8269342411124403  2.4198846981424063 
H  5.021411511262476  1.2739604515766971  1.1731555599388854 
H  1.0907497295862636  -4.502473766276894  0.48206976578556676 
H  0.416606860742472  -6.3166993192525895  -0.7372136488958377 
H  -0.22651901932495133  -6.071008580091195  -2.943912090703856 
H  -0.18872870101411188  -3.968836073120254  -3.922891816120807 
H  0.3924030815350402  -2.1080104735231386  -2.680436360502686 
H  3.4784812017351863  -1.9088707225571246  1.8351156772450639 
H  3.694020301824658  -2.471446025179486  4.062120251827682 
H  1.8411962541050215  -3.1640438980120957  5.2749674865359895 
H  -0.23545740622428224  -3.20138781469953  4.275733771668721 
H  -0.4795445156250562  -2.6662026477459553  2.055747147813237 
H  4.787046569049599  -3.4162083083885273  -1.613655285258185 
H  3.3382643705779835  -3.473464081572612  -2.20949034450603 
H  3.557509922251206  -3.635339524968754  -0.6647653790520331 
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