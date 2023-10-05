%chk=OJIKEZAW_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
C  2.7981711904012267  1.2830775321402221  0.1551400313592412 
C  1.1128724184777128  2.4042251371981704  -0.9602377508410371 
H  0.19997575491665653  2.565531507484808  -1.040688459750658 
C  3.244171677891522  2.800656394612986  -1.6647822667992214 
C  4.13456770585721  3.4775966282503132  -2.511079858350815 
H  3.809577208029303  4.107518882805281  -3.1151237139920855 
C  5.460251012275208  3.2196446941634136  -2.457284776926185 
H  6.041408471426738  3.684720445038291  -3.0140119305847346 
C  5.971034475769883  2.2581363382982174  -1.5686548511292904 
H  6.8776170350540795  2.060960141298535  -1.577847812485647 
C  5.152691458872109  1.6107307068489796  -0.6975151407139641 
H  5.504718211995185  0.9965414883450334  -0.09466814606343893 
C  3.749136507572705  1.8780170036670842  -0.7112218622835589 
C  3.2641449762890646  0.30406617876300923  1.1879294100923452 
C  2.8937797570946495  -1.0345502062718652  1.1006941641923254 
C  3.5495544659291363  -1.9871496925298437  1.9325956066926224 
H  3.3355229722454567  -2.888458211350163  1.8546671548059361 
C  4.486654555101529  -1.5968270386653405  2.8393154690279125 
H  4.958580731391496  -2.2410129912204733  3.315190087664951 
C  4.750993224893557  -0.24562342098909684  3.0688015932951966 
C  5.594675075714961  0.211946342799336  4.094784982956648 
H  6.025613256286288  -0.4175855540031441  4.626665162540499 
C  5.807504895601224  1.5227133305819711  4.345361800532439 
H  6.373709734531357  1.781236267073611  5.037095645070845 
C  5.1656194349847695  2.493971219178092  3.5483237027673455 
H  5.281323087895709  3.3970384041574144  3.7337939271069356 
C  4.36746976671925  2.1015896031445056  2.4953089286798456 
H  3.9765106883741335  2.74814955690385  1.9547796646669893 
C  4.138407471737624  0.7431933302595022  2.2242797774250667 
C  2.1886270896118556  -1.346799101358986  -1.6530681263946787 
C  1.4706043879276778  -0.6528213562529244  -2.631801828258318 
H  0.6025444482429506  -0.3707524724329776  -2.4550145665703633 
C  2.050622290873681  -0.3830127526860967  -3.860665852348638 
H  1.5732809488803148  0.07759476915212615  -4.513530006710057 
C  3.3348014153507526  -0.8005457857065166  -4.110592626812894 
H  3.729493612508152  -0.591066048215477  -4.925498234281513 
C  4.053962707228645  -1.5211532392396063  -3.1849168389815805 
H  4.913483270648284  -1.8251541384400851  -3.3723903420270496 
C  3.4470018336706576  -1.783017199777353  -1.9454949277095512 
H  3.9172078537316666  -2.266823848717438  -1.3065547699954254 
C  1.2238565242482882  -3.2438729833018067  0.39215261050844835 
C  1.5393321238251194  -4.2758368427055  -0.4764743744259917 
H  1.8926891485823643  -4.091670384780849  -1.3170496410322694 
C  1.3210432282453293  -5.580816810633287  -0.08098248908097457 
H  1.5122611438793532  -6.278236572680812  -0.6667095537161684 
C  0.8205333516675389  -5.863815634227369  1.1817343316244784 
H  0.7016899525783641  -6.748172373622641  1.4465462535220843 
C  0.49704360031214034  -4.8410771534254655  2.044282191545925 
H  0.15271428770357032  -5.023789981820396  2.8881540282714058 
C  0.6922007203649629  -3.550011504205642  1.636047751674304 
H  0.4627058916846809  -2.8555370497814057  2.2097909483757316 
N  1.5007227975212478  1.4966399316468877  -1.9828502194830842e-16 
N  1.887121755544733  3.050577018198384  -1.7592610483604387 
P  1.5007227975212478  -1.4966399316468877  2.220446049250313e-16 
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

-Pd 0
lanl2dz


