%chk=SOQURAWA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
C  2.653976650971198  1.2830775321402232  0.1551400313592399 
C  0.9686778790476847  2.4042251371981704  -0.9602377508410348 
H  0.055781215486629954  2.565531507484809  -1.0406884597506585 
C  3.0999771384614867  2.8006563946129877  -1.6647822667992118 
C  3.9903731664271573  3.477596628250314  -2.5110798583508007 
H  3.665382668599271  4.107518882805282  -3.115123713992073 
C  5.316056472845174  3.2196446941634145  -2.457284776926178 
H  5.897213931996713  3.684720445038292  -3.0140119305847364 
C  5.82683993633985  2.2581363382982182  -1.568654851129285 
H  6.733422495624054  2.060960141298536  -1.5778478124856559 
C  5.0084969194420745  1.6107307068489802  -0.6975151407139593 
H  5.360523672565153  0.9965414883450341  -0.09466814606344075 
C  3.6049419681426604  1.878017003667085  -0.7112218622835513 
C  3.1199504368590265  0.30406617876300923  1.1879294100923372 
C  2.749585217664622  -1.034550206271865  1.100694164192325 
C  3.405359926499097  -1.987149692529844  1.9325956066926153 
H  3.191328432815422  -2.888458211350163  1.8546671548059288 
C  4.342460015671491  -1.5968270386653403  2.8393154690279037 
H  4.81438619196147  -2.2410129912204737  3.315190087664952 
C  4.606798685463526  -0.2456234209890964  3.068801593295195 
C  5.4504805362849345  0.21194634279933644  4.094784982956643 
H  5.881418716856255  -0.4175855540031441  4.626665162540502 
C  5.6633103561711975  1.522713330581971  4.345361800532434 
H  6.229515195101317  1.7812362670736108  5.037095645070827 
C  5.021424895554728  2.493971219178092  3.548323702767335 
H  5.137128548465682  3.3970384041574153  3.7337939271069356 
C  4.223275227289213  2.1015896031445056  2.495308928679842 
H  3.832316148944085  2.7481495569038508  1.9547796646669748 
C  3.9942129323075877  0.743193330259502  2.224279777425055 
C  2.0444325501818286  -1.346799101358986  -1.653068126394679 
C  1.3264098484976523  -0.6528213562529244  -2.631801828258303 
H  0.45834990881292315  -0.3707524724329776  -2.4550145665703527 
C  1.9064277514436558  -0.38301275268609714  -3.86066585234863 
H  1.4290864094502904  0.07759476915212571  -4.51353000671005 
C  3.1906068759207216  -0.8005457857065172  -4.1105926268128945 
H  3.585299073078125  -0.5910660482154776  -4.925498234281513 
C  3.9097681677986156  -1.5211532392396072  -3.1849168389815787 
H  4.7692887312182535  -1.8251541384400858  -3.3723903420270487 
C  3.3028072942406186  -1.7830171997773538  -1.9454949277095452 
H  3.7730133143016227  -2.266823848717439  -1.3065547699954165 
C  1.0796619848182607  -3.243872983301807  0.3921526105084483 
C  1.3951375843950913  -4.275836842705501  -0.47647437442599244 
H  1.7484946091523363  -4.0916703847808495  -1.3170496410322694 
C  1.1768486888153014  -5.580816810633288  -0.0809824890809755 
H  1.3680666044493253  -6.278236572680813  -0.6667095537161674 
C  0.6763388122375119  -5.86381563422737  1.1817343316244784 
H  0.5574954131483372  -6.748172373622642  1.446546253522086 
C  0.3528490608821193  -4.841077153425466  2.0442821915459097 
H  0.008519748273549066  -5.023789981820397  2.888154028271388 
C  0.5480061809349344  -3.550011504205643  1.6360477516743037 
H  0.31851135225465965  -2.8555370497814057  2.2097909483757228 
N  1.3565282580912201  1.496639931646888  -1.9828502194830842e-16 
N  1.7429272161147042  3.0505770181983856  -1.759261048360438 
P  1.3565282580912201  -1.496639931646888  2.220446049250313e-16 
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


