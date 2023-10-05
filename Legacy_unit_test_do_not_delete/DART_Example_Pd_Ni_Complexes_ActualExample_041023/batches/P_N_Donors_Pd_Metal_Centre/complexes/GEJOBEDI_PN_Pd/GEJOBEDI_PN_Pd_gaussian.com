%chk=GEJOBEDI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
C  3.3941038108942374  0.9539552521268835  2.003147485687358 
C  2.598604029095547  1.6478390128220903  2.899839296624495 
H  1.8197534075322614  2.0608156774732245  2.60324886565257 
C  2.9550488145720646  1.7295463653611245  4.238261126036922 
H  2.415177585046596  2.195976207417092  4.8353497191036015 
C  4.088650727676248  1.1302928483549621  4.676901496912999 
H  4.337234816530952  1.198241706575267  5.570046632618764 
C  4.872059062290421  0.4188945503844419  3.7947261555540868 
H  5.635855348357253  -0.010203755821134308  4.104359579140852 
C  4.5373432650491665  0.3335601654587015  2.459851960048554 
H  5.080105329949427  -0.1376708915749405  1.871414459720628 
C  2.9922238016773255  -0.6327086287124573  -0.48074051449134036 
H  3.773373208432309  -1.1154367669902385  -0.17112481908624932 
H  3.0436096569160265  -0.5569438508043095  -1.4457220854604094 
C  4.211275121048454  1.999272551866862  -0.5300440216261441 
C  4.2223169833237035  3.3703839900680275  -0.309796635009572 
H  3.5799775987508182  3.757936630223362  0.23871568585034952 
C  5.190288262567515  4.158446985388203  -0.906793559596409 
H  5.202525336072781  5.074231357358906  -0.7494547910546137 
C  6.132246030631001  3.5945037278374676  -1.7309596929991404 
H  6.770688358002478  4.134902048008168  -2.1379051809667518 
C  6.14133775097687  2.2349233278911043  -1.9600718417447487 
H  6.779837588153289  1.860074865755875  -2.523310049220578 
C  5.186730804154515  1.4213830262013918  -1.341835532066031 
H  5.202938335272266  0.5005655865197491  -1.4712395069401252 
C  1.9371338475342332  -2.4283840954407028  1.5211808543083254 
C  3.239674369173134  -2.885048833255502  1.7161289672132556 
H  3.8728003928354973  -2.7764787152765695  1.0446051719211393 
C  3.5948807695434253  -3.4992964365676116  2.9066262015491224 
H  4.474561963681396  -3.7650807746541983  3.041628654519063 
C  2.6587957720702793  -3.7190941051633843  3.8890236885092744 
H  2.9044665103937266  -4.118980857099203  4.693016234468395 
C  1.3508113413019547  -3.3397440318635603  3.6717192229427513 
H  0.700286102490198  -3.5294785551999066  4.308989600665794 
C  1.0031770777309432  -2.6786256970868747  2.510031835469872 
H  0.12558003319653288  -2.397341997589807  2.390166294876558 
C  1.3679772172197635  -2.8550957908331425  -1.258976195065787 
C  1.9574482212603752  -2.751423016805064  -2.493845439690777 
H  2.4486637419974735  -1.9910433538509733  -2.709401674471425 
C  1.826790233384918  -3.7708605316107353  -3.422937078664655 
H  2.241010997673328  -3.694801070723106  -4.252869947541332 
C  1.098860696748659  -4.878475754172207  -3.1330435664181424 
H  1.0132943514460473  -5.556102866495532  -3.7639529079212917 
C  0.4902898169057843  -4.999456475194032  -1.910028844778083 
H  -0.014215889411786309  -5.755972097558777  -1.7116228332484589 
C  0.6289354458869304  -3.9849599335579025  -0.9592733101779585 
H  0.22530988097250448  -4.0692642022580685  -0.1261008588652363 
N  1.4622056387526345  1.5342928892489853  0.0 
P  2.9427023051627987  0.9982861425824154  0.2616818837374227 
P  1.4622056387526345  -1.5342928892489858  0.0 
H  1.3568674991382481  2.241176755816913  0.08187481620472076 
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

-Pd 0
lanl2dz


