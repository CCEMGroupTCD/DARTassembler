%chk=SUWODUZI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  2.286703494418004  -0.9302449660989522  -1.242914027636063 
P  1.4713292748395923  1.5255458580455712  2.429184246265767e-16 
N  1.471329274839592  -1.5255458580455716  -3.3306690738754696e-16 
C  2.9064457915649995  0.7051175448623357  -0.7945254095051314 
H  3.199182897026751  1.204181896802162  -1.5976680927867482 
H  3.666514470132875  0.6337913999115692  -0.1642662922314574 
C  3.708031741245269  -1.9192383267658784  -1.750869465964806 
C  4.970738896758723  -1.3246080308518886  -1.8360719342179954 
H  5.074900034452462  -0.3979719729814015  -1.6510746270792285 
C  6.0731240963920055  -2.090782439071536  -2.1910110104251608 
H  6.93063344422289  -1.688129790670028  -2.255426715946929 
C  5.922233017070039  -3.441111997593704  -2.451140335112674 
H  6.681158143455244  -3.972641830488297  -2.663776595624003 
C  4.670231493787996  -4.015950885834979  -2.402912041912619 
H  4.568738141236686  -4.9388483081410754  -2.605547249338819 
C  3.563815063602785  -3.265111781943694  -2.06154417130384 
H  2.7043039051904127  -3.668442971075092  -2.0401493018603185 
C  1.1542602340925465  -0.7514652863131133  -2.6288116025813593 
C  0.03366896175396783  -1.5724628711412616  -2.700810394689157 
H  -0.1346154510699673  -2.2141362465021994  -2.0193772169281914 
C  -0.8403081233024603  -1.4490044519749636  -3.779245851048632 
H  -1.6048138842213955  -2.008708793536942  -3.834168356075441 
C  -0.5981610625617468  -0.5157539485623881  -4.765813405174981 
H  -1.2013221436024184  -0.4304547624948445  -5.49567985129646 
C  0.5192672428821379  0.29589466787862073  -4.692442596114164 
H  0.6833136567446187  0.9358290689661799  -5.377062517240896 
C  1.4023157809757758  0.1841643556752437  -3.631143657687457 
H  2.171197145415751  0.7398904228589118  -3.5865438183434013 
C  0.9273606256415036  2.704572667049086  -1.2649555535040093 
C  -0.26646616134898227  2.5191073278671627  -1.9573407285519724 
H  -0.8503860798971699  1.8089071432669854  -1.715243461907813 
C  -0.602085922841505  3.3651469720988496  -2.9956342448969555 
H  -1.411834401493343  3.2315259774070535  -3.4728120850138424 
C  0.22334726597750665  4.389891204961694  -3.3364614562672363 
H  -0.018136769691922794  4.969391578770577  -4.049765569021527 
C  1.4104480457589113  4.598421432566623  -2.662196484006724 
H  1.9782558890140298  5.318376718871225  -2.911357674524117 
C  1.7742934656977516  3.7611508233196576  -1.6230806351810032 
H  2.59048701125138  3.900718823768215  -1.1571826104509046 
C  2.192485228915586  2.5185061823852593  1.3441212126099384 
C  1.3389244805264124  3.1481945230757704  2.2557659254272755 
H  0.3974919392286276  3.0677844492959148  2.1545563892148 
C  1.8502852314107148  3.880401463371257  3.2950509436192164 
H  1.2610071040779582  4.298897499811564  3.911026573113568 
C  3.211687488865774  4.013664543551917  3.451479473168876 
H  3.5616770817857555  4.518697873711555  4.175345262024532 
C  4.0656055969091875  3.406952208345044  2.5483421375086723 
H  5.005970118365075  3.5040083926735113  2.6498620560994737 
C  3.562456720770351  2.658854352760353  1.497337310255329 
H  4.1553035174655015  2.243476374013896  0.8828691008084651 
C  2.210563330264061  -2.269443492466011  1.0709979742298277 
H  2.947026405299243  -2.769250054149339  0.6147299827949443 
C  1.3214167524869203  -3.303624645840995  1.7145191045986885 
H  0.5797119232548655  -2.856194706976705  2.1724584413036463 
H  1.8409027808550205  -3.822271889925587  2.3635274269905038 
H  0.9651807425847894  -3.903352104163437  1.0258280846317285 
C  2.873078837647996  -1.3595688153597556  2.1007769478142913 
C  2.1316508738490434  -0.6255843670426662  3.012278294014967 
H  1.182493865653335  -0.6720107721398799  2.988341863563879 
C  2.75064930818031  0.17371890107552712  3.9559273780875603 
H  2.226400321464773  0.6903198513803965  4.556622407259473 
C  4.123016985522356  0.22310645609568547  4.026856692168612 
H  4.545846993497204  0.7514664170559153  4.693494146869505 
C  4.8918540273798525  -0.4983542684019988  3.12332019258759 
H  5.840000874795925  -0.4592207938590025  3.1651379747529877 
C  4.2678789056535145  -1.2771971453524116  2.158113395838719 
H  4.79402138889239  -1.7594129903368794  1.5300457671734993 
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

-Pd 0
lanl2dz

