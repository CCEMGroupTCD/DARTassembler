%chk=IVIJIMOT_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.6231655591774976  1.3628769451054636  1.2982597836612787e-15 
N  1.623165559177498  -1.3628769451054636  0.0 
N  3.8889881679952527  -1.7064837884685402  0.45895434692836534 
C  2.8291613234673685  -0.9083800352967306  0.2901753530009741 
C  1.7178462601766091  3.133518468293944  0.29511567430263425 
C  2.24402952554874  3.9748711040678715  -0.7203669959180894 
C  2.968303831290827  0.5506260893142696  0.4969196174055588 
C  2.5168328447431874  3.6421905736441116  -2.199405204928744 
C  2.2048509227432147  5.81256903343316  0.8782633983534512 
C  1.6417748262939083  4.962272226254926  1.8069914098347177 
H  1.4445875698914081  5.308224134600925  2.6486635924712862 
C  1.3521992099082105  3.6444310180297657  1.5705438420770879 
C  -0.7775660350843172  2.4374971301924364  2.091143008050114 
H  -1.2709815691113395  3.20891983550399  1.8033451961914269 
H  -1.2739641309519498  1.9760855792388479  2.7712113289309075 
H  -0.6355251838103066  1.8502664210119621  1.346307307533043 
C  5.316553797797646  -1.3314377787494123  0.4168625044800778 
H  5.812607631711261  -2.0667448793640255  0.024840934097357437 
H  5.629988155993257  -1.222415456129888  1.3279168037801719 
C  1.4589703801444698  -2.7927209798133847  -0.18530607003681682 
H  1.7285685148426435  -3.038817981603994  -1.0847938417585616 
H  0.5268050589082376  -3.0332238155075006  -0.07106072042616635 
C  0.5739206528483194  2.879319297015725  2.657498722459761 
C  2.3054103325823476  -3.53115882529482  0.8223959917243658 
H  2.0344530671660017  -3.291955182463535  1.7226802343801846 
H  2.191570595062728  -4.487512580988545  0.7126458293920359 
C  2.4634862733623466  7.29608882754307  1.2290027686759364 
C  2.4664607218840358  5.304084212049396  -0.36209132486723145 
H  2.81381983000744  5.877950816667527  -1.0052096655385174 
C  3.7434179266476377  -3.154245952239911  0.5995831206058708 
H  4.277105709102558  -3.4612616197238313  1.3484896172913894 
H  4.071655906250048  -3.589745526967891  -0.20111976383896712 
C  5.643551584319158  -0.09241802273664179  -0.34356615872252005 
H  5.1056042676047  -0.06951732351204365  -1.1497198490414073 
H  6.575357890134655  -0.1272561349891954  -0.6121321922733238 
C  5.411390501894026  1.184041004798952  0.43317658043253116 
H  6.159470509347237  1.3341665577586626  1.0301950314014994 
H  5.370386209950256  1.9308589318725065  -0.18598743809871077 
C  4.143424346236317  1.139524050999205  1.2324485672624406 
H  3.916529659132207  2.042037113104617  1.5077467032769083 
H  4.299700340895826  0.6190322651871649  2.035444335866103 
C  0.24614503218611206  3.784858542090048  3.8609540054386313 
H  1.0619262451676168  4.10934993992918  4.248524947686446 
H  -0.24131123692954803  3.280817992085114  4.517538881635609 
H  -0.28673759478358707  4.525988842685042  3.5668064856227253 
C  1.313492352778758  1.6721428194118069  3.2303681251876446 
H  1.4396014299405442  1.0181783597437106  2.5394288542322174 
H  0.7975691876040856  1.2906757604325612  3.942388595785173 
H  2.1680209834216653  1.9508218126723542  3.566790253066986 
C  2.562839420661653  2.227542650422542  -2.619065196227856 
H  3.2850488379564338  1.7849956529360385  -2.1681489372767984 
H  2.699765321342335  2.1809998534166897  -3.568046479378203 
H  1.7344355530319646  1.7998263719327836  -2.3930820393089265 
C  3.2617635034089663  7.371876594275015  2.5227039468476344 
H  3.9664592059043056  6.719333084610479  2.5006779422326337 
H  2.682820161115624  7.192450589690958  3.267228872752221 
H  3.6404453162783392  8.248806271816457  2.6148217260709754 
C  3.817661176991346  4.245968924987727  -2.6334993219924216 
H  4.53382220523515  3.869479879632062  -2.1169929422138423 
H  3.7902005055066113  5.197080718062846  -2.497986080591306 
H  3.9622876850930435  4.059635367206141  -3.563976177372932 
C  1.1647062272390163  7.979070286979212  1.537896031565226 
H  1.3389951429404872  8.817903289551968  1.972762387654519 
H  0.6397575630261484  7.422510923550181  2.1155224722001753 
H  0.684105238812426  8.136199656412805  0.7209793153787991 
C  1.4309293126129268  4.2787444519055615  -2.994149169345813 
H  1.598545523742969  4.145662114097018  -3.9304946809652503 
H  1.4064854058127063  5.220132252622022  -2.8054397522978842 
H  0.5884028268688424  3.8832045449932147  -2.7594917533975787 
C  3.1812740590015522  8.013248350646139  0.18942365518347404 
H  2.597483028839016  8.164591427319902  -0.5584833218161458 
H  3.9361251275169904  7.494580584976877  -0.09485767240725806 
H  3.483829482323297  8.856681005463738  0.532934517611335 
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


