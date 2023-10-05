%chk=EDONIDOQ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5049920639325645  1.4923467718663783  1.3672291662134427e-15 
N  1.5049920639325642  -1.4923467718663783  1.1102230246251565e-16 
N  2.3113986781839544  -1.7521237267326777  -1.0339292752466376 
C  1.7510166863280383  -2.398039845674977  0.961462079975243 
C  2.7297640273565227  -3.2495486799191546  0.5073310772526582 
H  3.091067293412465  -3.9733190799108202  0.9657052672254622 
C  3.0753404230822983  -2.827365118512682  -0.7587411581254735 
C  1.0315441046398512  -2.394678756888805  2.2366646392892995 
H  1.0319783341855715  -1.5075782049545277  2.60156217586993 
H  1.4642420488092258  -2.9929707234961107  2.8500422895218196 
H  0.1267125049401585  -2.681055451003613  2.0953428247920773 
C  4.0214445470609945  -3.3767147639538058  -1.724181582825124 
H  4.019240358580236  -2.8361744115321317  -2.5178711499463167 
H  3.765781128893668  -4.2740781450701  -1.9482963945679872 
H  4.902303479552545  -3.3799589797696408  -1.3428628105182754 
C  2.4888045265234817  -0.7418687890585853  -2.062276578388951 
H  1.6284405301606795  -0.34636131644762314  -2.272004718779149 
H  2.8277826412850677  -1.163230718373165  -2.867276598411942 
C  3.448359966232421  0.35402715874089496  -1.6299577033422104 
C  3.1995552488164383  0.8894103451840651  -0.24769521153194132 
H  3.818673057534542  1.615768920439561  -0.07788069777970708 
H  3.3807608482796914  0.18884804310431047  0.3977016366381542 
C  1.37327368848818  2.896042206125384  -1.139276198913476 
C  2.0096448777167826  4.08659522871676  -0.8079296885677051 
H  2.45333370092908  4.165810195525148  0.00476255157615499 
C  1.9834794617544853  5.158326166292177  -1.6926273422614706 
H  2.4091622135882718  5.9549379809264  -1.4701367321755785 
C  1.3218453744370233  5.039418960309432  -2.9090153317060468 
H  1.3037584286571602  5.756081655711275  -3.501632428499341 
C  0.6853742086534992  3.8488637766921725  -3.240361977815873 
H  0.24177745872494683  3.770068249930697  -4.053876162335894 
C  0.711537869552568  2.7772100430228264  -2.3556007908502004 
H  0.28596416392924784  1.980178215418551  -2.578050193789002 
C  1.482974698072097  2.164239564414852  1.614879998923156 
C  2.429985671079401  1.8041832724526308  2.5665044698265715 
H  3.1203942731176824  1.2240255495145456  2.3386262814301007 
C  2.344663485727834  2.3108183181762394  3.858206395632794 
H  2.978143809886685  2.070309307960657  4.495825580498501 
C  1.3123320824321845  3.1774324519559975  4.19822031726369 
H  1.2554272998510345  3.5160863020715043  5.062251506238617 
C  0.36542381485888886  3.5373863347443706  3.2462509348427115 
H  -0.3258673198809616  4.116788723854387  3.4731338309007924 
C  0.450745026394588  3.030778655314423  1.9548305230461835 
H  -0.1836456171379104  3.2717314397334833  1.3175913159354438 
H  4.387414471638605  0.07664857435041417  -1.733166456836856 
H  3.3815332453996803  1.0975191010408487  -2.334120512014891 
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


