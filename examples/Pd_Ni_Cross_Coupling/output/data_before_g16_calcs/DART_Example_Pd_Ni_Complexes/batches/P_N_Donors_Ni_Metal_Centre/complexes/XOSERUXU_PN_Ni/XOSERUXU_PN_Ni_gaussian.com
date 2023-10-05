%chk=XOSERUXU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4274310000837167  -1.4291748458463722  2.614285692557143e-16 
N  1.4274310000837165  1.4291748458463718  -3.1380131811691166e-16 
C  1.4078318186803611  -0.4512983641396479  -1.5515587222574736 
H  2.070755839447268  -0.8207385145414812  -2.1877175530232944 
H  0.5120182188358807  -0.5046573847112972  -1.9692974315690783 
C  1.740866179578281  0.9780961749100527  -1.243095262190945 
C  2.329125141413255  1.8161877656493666  -2.183644883079135 
H  2.546103971355798  1.4857068161891274  -3.046893844503734 
C  2.5965379633105354  3.132306556479381  -1.856245782425097 
H  3.0065908463693627  3.7113326094532773  -2.4883748776838246 
C  2.258815603232107  3.5963997511837746  -0.6016194738645863 
H  2.4260977778490442  4.499256482639784  -0.35749488165806687 
C  1.6751460964013172  2.7233284166116523  0.2938838994090046 
H  1.4373862016126784  3.04817359657636  1.154163779148038 
C  3.1701530400567246  -1.8599332878866948  0.3046549950409566 
C  3.451500534868592  -2.7528190000224178  1.3413787539213013 
H  2.7478967508112286  -3.1090108553993074  1.8709165780881294 
C  4.771923668507571  -3.117528945425439  1.5940671269723976 
H  4.965808376356936  -3.7362569845943083  2.288830067556027 
C  5.805276601656246  -2.5836014471786455  0.8389021104568153 
H  6.705256397093961  -2.8312303716790708  1.020261569969357 
C  5.526018164659732  -1.6962564706800016  -0.1724321607161307 
H  6.235623189265609  -1.3308784124825426  -0.6878455808382509 
C  4.209874126297192  -1.3259534937292794  -0.4478767931637398 
H  4.025276087224341  -0.7107646856171962  -1.1473046696056965 
C  0.6846367319260431  -3.0304742786283048  -0.4275933719836988 
C  -0.15087330490362194  -3.6734149675653778  0.4768929512219391 
H  -0.37820700204518776  -3.2488036670210962  1.2951206540003624 
C  -0.6556545472370656  -4.939398290220701  0.18930299286647448 
H  -1.2235704503903961  -5.38025308264246  0.8110114526322446 
C  -0.3278283508907105  -5.549996806235869  -1.0023162782219035 
H  -0.670744369000376  -6.413956393537989  -1.2001292606890954 
C  0.4936853504216622  -4.917983289747561  -1.9071017575108662 
H  0.7081563753158874  -5.343859341071115  -2.7286039172097287 
C  1.0135976765562449  -3.6573992014983197  -1.6267366569969912 
H  1.5889108159122314  -3.2284903512050693  -2.2489091891972666 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz


