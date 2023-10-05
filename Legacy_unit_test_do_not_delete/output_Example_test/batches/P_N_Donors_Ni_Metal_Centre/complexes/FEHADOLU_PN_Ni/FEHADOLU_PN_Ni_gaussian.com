%chk=FEHADOLU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4620502128860016  -1.3937392779856648  -2.113811983227793e-15 
C  4.2188061568055435  -1.0763382102340593  0.45987263220882846 
C  2.8521048051542843  -0.5132406155861793  0.21910715649339224 
C  2.655840672079349  0.944075004402319  0.16058695651181443 
C  3.8556805494091844  1.7990228406447994  0.31752460678572236 
C  4.579683216493427  2.25419833151342  -0.798092276632372 
C  5.78037557582453  2.9654729799776205  -0.6306029256394191 
C  6.232088520928514  3.2679146716613166  0.600075914821026 
C  5.533330685901596  2.8375798131454273  1.7378192855046892 
C  4.374763669658711  2.100611962541044  1.5663897779815295 
N  1.4620502128860002  1.3937392779856643  -1.7068383456310868e-16 
C  1.69495988614825  -3.1677228085161384  0.1833448322356377 
C  2.0583205004098004  -3.9530114290548624  -0.9592925627658984 
C  2.334542627902362  -5.30790678848623  -0.7455999768513779 
C  2.2084971180486432  -5.931711282434735  0.5109350297343584 
C  1.8627088854700675  -5.152284009229067  1.5625682709327282 
C  1.5829755805352144  -3.77972701078902  1.455043138319567 
C  2.1584281142634243  -3.4211485966667627  -2.416077753922345 
C  0.917575379146924  -2.6253312135167075  -2.8319353646861396 
C  3.4179391177863527  -2.5269490755043695  -2.6080614729134943 
C  2.2709760656797373  -4.549861756902592  -3.4590397573253817 
C  2.501706628004721  -7.425940485051369  0.6667877593189907 
C  2.905288152414368  -8.106741388769432  -0.5412034834687597 
C  3.4927771621709915  -7.633610068288137  1.7344364236771421 
C  1.2125811239021653  -8.048640267362394  1.2969168459212357 
C  1.2299546380872406  -3.040549783546812  2.78693483471223 
C  1.1465155658826172  -4.042206038092261  3.984267235910683 
C  -0.1878615629072009  -2.419413776494491  2.641428714769322 
C  2.2424026960065886  -2.0169256003624763  3.1820584123630242 
C  0.8592811600234663  3.472459079286837  1.1401574437724173 
C  0.40866101792953136  4.761912927209811  1.0318479261117468 
C  0.2456618909335906  5.394944638330946  -0.2118579584580928 
C  0.5887168152091418  4.735515059559002  -1.3740874579851745 
C  1.0537715417132336  3.4234831904113707  -1.3298227350426268 
C  1.175185141746036  2.806525690840318  -0.0740672168583908 
C  0.9794684017236779  2.8023117068529477  2.477358302339492 
C  1.6606356721004367  3.6645093603027075  3.5044536319721535 
C  -0.40462057169983257  2.3322541068785445  2.9841559506705324 
C  1.4173330566136093  2.681342403150986  -2.603125749555648 
C  2.1824251917943296  3.565762964779644  -3.6061665038887356 
C  0.18312816207356097  2.1036132376475627  -3.277546263493477 
H  4.165627389563278  -2.0330175376094597  0.4814914501031029 
H  4.8077710342969295  -0.7925966695841105  -0.24308706165077423 
H  4.55440647483816  -0.765177911568457  1.2998707542726626 
H  4.244784539996245  2.068930617473538  -1.6840649007099475 
H  6.28779883862089  3.244396987602523  -1.3985424351996272 
H  7.0436029845653065  3.784211425556347  0.6951097889173614 
H  5.852816539961443  3.052758878367327  2.615914348095242 
H  3.9184588175622617  1.7562031820883928  2.3434923991972685 
H  2.6117629568908254  -5.83592471237197  -1.5078177352714575 
H  1.7925647275529362  -5.580562328157338  2.4286746492986393 
H  0.8044945648706296  -1.897354571094468  -2.2241978658181947 
H  1.0192682875695405  -2.2762107842616417  -3.724235727934586 
H  0.15495009203310728  -3.1982793521055193  -2.789853160107572 
H  3.0480170890789937  -5.086804194286994  -3.2410914907795028 
H  1.4913971788898226  -5.105839268076652  -3.413012522560275 
H  2.3557067101961096  -4.183863364622765  -4.347431672376613 
H  3.374903993459375  -1.785240456591046  -1.9929257214700682 
H  4.191562895262599  -3.067539698801527  -2.417512899246327 
H  3.4578346679248533  -2.1882661077098327  -3.5073978796073577 
H  2.2266909372756842  -7.966856169145158  -1.2112833928353288 
H  3.741640213075182  -7.761575956756005  -0.8588704870612545 
H  2.993229578206477  -9.048527066143087  -0.36576057673413176 
H  4.338739449101675  -7.28678169971505  1.4370701144826357 
H  3.2177758415854134  -7.165068637825044  2.5183076513770755 
H  3.5799510481147827  -8.576257825836313  1.9387605994343422 
H  0.9910149101448419  -7.573797966902518  2.100282903590612 
H  0.48377495507490376  -7.976427912735021  0.6710085825802541 
H  1.3536449974370734  -8.977781485131333  1.5026490772953773 
H  -0.802122496784327  -3.105412409884398  2.394032190744405 
H  -0.46522748182933005  -2.04499291224633  3.4822162642942573 
H  -0.1855526110013519  -1.7447349431914436  1.95943442547215 
H  2.3218809237788136  -1.360256893053925  2.4795700208917446 
H  1.9730191514682527  -1.5860852204688585  3.990665067656631 
H  3.08848658735456  -2.442002700045171  3.324993356050222 
H  2.003293593674984  -4.450811294466372  4.120395284290475 
H  0.878493933627754  -3.5691690752875287  4.776693236534047 
H  0.5066320157940762  -4.728540928777492  3.786708332911611 
H  0.19106246419163408  5.2412449016702976  1.839029806279418 
H  -0.09703091882248183  6.286582306601341  -0.25870864322176895 
H  0.5105950552893778  5.167915630827252  -2.226180819415706 
H  1.5210249896658032  2.0139915114230393  2.3654737120965748 
H  2.520752207351398  3.922316925659775  3.1829339824442227 
H  1.1198600199777564  4.437691078735331  3.65449993933378 
H  1.765729732687056  3.1803655891841744  4.331153563819584 
H  -0.8139870384507739  1.7622618080691443  2.324947239261275 
H  -0.3150076593125595  1.8276996101319627  3.8074520126957783 
H  -0.9507852256262972  3.1014346824230525  3.140419935665133 
H  1.9902759886298529  1.939600951698088  -2.374901154366613 
H  2.9731767014071853  3.8962075409457295  -3.174579754639604 
H  2.431502465678468  3.073660172791705  -4.388454682089137 
H  1.62110190506417  4.29476194690523  -3.863911148516559 
H  -0.28049146469346575  1.5278014357731988  -2.6675696897641115 
H  -0.38336440536558425  2.8292486530903083  -3.531013647768933 
H  0.42703615524872696  1.6081468789767817  -4.055557181341559 
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