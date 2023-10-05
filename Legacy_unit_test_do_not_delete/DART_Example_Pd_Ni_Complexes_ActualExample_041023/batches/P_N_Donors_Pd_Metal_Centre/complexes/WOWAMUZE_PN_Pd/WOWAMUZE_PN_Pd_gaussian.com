%chk=WOWAMUZE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.454426734146482  -1.5416688603587994  1.5611738195316817e-15 
N  1.4544267341464818  1.5416688603587994  -1.8879998351835515e-16 
C  2.6732242030111966  1.4041092323134854  0.37588882481531644 
H  3.2334372823438384  2.1251578657672923  0.20369131049222847 
C  3.3025336174567967  0.24948013149236595  1.046452625230423 
C  2.9191936669645644  -1.0392495731068456  0.9687073305803037 
C  3.717533514397279  -2.14541072343527  1.6180463695308391 
H  3.732568285354878  -2.9271886531779114  1.0449704284640713 
H  3.314654259076173  -2.3933238809411224  2.4651956963719868 
C  5.219554343671896  -1.6045225038950954  1.8555489008603396 
H  5.617070651928653  -2.180839210320117  2.525780419490521 
H  5.7035218741409945  -1.7677785969976207  1.0306186871352312 
C  5.449377512010894  -0.6112022557040564  2.1457696544744147 
H  6.327401011119737  -0.41854279579197473  1.7828228849350454 
H  5.548525807122936  -0.6282655276403372  3.110423955417091 
C  4.542534809058335  0.6732003069442914  1.8130349415812324 
H  4.281692439712488  1.1128074187090504  2.637438887253435 
H  5.052383297807883  1.3059659027759178  1.2824368646419653 
C  1.1132277951841911  2.800788026551355  -0.6340951816388322 
C  0.5413177841224011  3.7998230979619216  0.15752858337773767 
C  0.1255408497281556  4.962910581341903  -0.48751276167741026 
H  -0.26426148576446096  5.644007194888095  0.009256986053274272 
C  0.2803316774325262  5.1207881912868505  -1.837686483809752 
H  -0.00940406564369467  5.903420205219589  -2.248085838038851 
C  0.8679984470797721  4.115810627791527  -2.60317231955289 
H  0.9791289129643856  4.239558956570509  -3.5182196520317692 
C  1.290673580087779  2.929158372537289  -2.010931129633879 
C  1.8997829584416732  1.8188571080999805  -2.8288741150492553 
H  2.1264789912773217  1.0839678500161656  -2.2530197089277118 
H  2.6918610763279016  2.141152017062702  -3.2648799141191214 
H  1.2675934712735017  1.52532354091443  -3.4892628393998324 
C  0.3780202728011548  3.6216949395348994  1.6367717892001465 
H  0.7194353708968446  2.7641503111168917  1.8952468384094212 
H  -0.5537555905526117  3.678861567379667  1.8644432128179709 
H  0.8632643451730926  4.31001758279892  2.097269942288663 
C  0.884259681404993  -3.089780443441913  0.7397209462304033 
C  -0.12951245513218956  -3.0262002560745764  1.6920839783537307 
H  -0.4955417953320367  -2.205036219132579  1.9266081701676943 
C  -0.5943641791412901  -4.193173625817085  2.2935354635107625 
H  -1.2669338606637341  -4.15170740590192  2.9358954841400076 
C  -0.056110870361964826  -5.416599987663062  1.9362933569101732 
H  -0.37380353883376505  -6.194931924471293  2.3328375972461304 
C  0.9476932168173774  -5.488401891980075  0.9938445248251124 
H  1.306291380663457  -6.312567364651223  0.7602471563832395 
C  1.4242264695320852  -4.332841884375467  0.39252362224314075 
H  2.101846234692542  -4.384312444649696  -0.24107455901131805 
C  2.093604608372159  -1.9649229045821077  -1.6391648633218732 
C  1.1876353048140427  -2.4226621591946462  -2.595628247773779 
H  0.3079288998055745  -2.601747728151228  -2.352202427494624 
C  1.5904561088831473  -2.6132716717534006  -3.9068918648741984 
H  0.9836609718131077  -2.917343889889077  -4.54274645342902 
C  2.910698096023304  -2.344004888416161  -4.268160289710188 
H  3.1847936583576395  -2.46202114967583  -5.14873342397386 
C  3.8100299690421284  -1.900447505515627  -3.322475542954209 
H  4.690718844232286  -1.729241851832922  -3.565005993131875 
C  3.4055132162743087  -1.7106833009432434  -2.012334072198761 
H  4.016800204790912  -1.4099633493888597  -1.3784883587893888 
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


