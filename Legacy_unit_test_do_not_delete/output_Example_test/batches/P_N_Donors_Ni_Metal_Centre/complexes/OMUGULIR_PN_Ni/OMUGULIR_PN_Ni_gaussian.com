%chk=OMUGULIR_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4468467213219236  1.4095157200258532  1.5672969722864733e-15 
C  2.9252708546447805  0.4312892374047792  0.4869661597148412 
H  3.7150747628327236  0.7736307119582875  0.04051108760615052 
H  3.063920544521581  0.4900118957079309  1.4463245926787032 
C  2.674929711894015  -1.0086132171508697  0.0808984638673774 
N  1.4468467213219236  -1.4095157200258537  -2.7755575615628914e-16 
C  3.9110838316205054  -1.8253386879238276  -0.17670411312654122 
H  4.6660480400105415  -1.226486686482852  -0.28887681221219835 
H  3.7984528037041034  -2.317247628809767  -1.003347310296256 
P  4.296989103352071  -3.0378829048614375  1.2074293534224665 
N  1.5900208953711181  2.825693111054444  0.8729115264588929 
N  1.5336039429546968  1.7449740858192675  -1.6108295562578336 
N  5.101411730322503  -2.072471749337015  2.336148511586741 
N  5.380117616644163  -4.032871304121204  0.3621453898739901 
C  1.1598214062159509  -2.8125461062096457  -0.4724785539400525 
H  1.9966104388844441  -3.1859926151924225  -0.8179006576542169 
C  0.716501347405603  -3.6778090668651187  0.6742061387170637 
H  1.392182112872166  -3.6729079260677153  1.35765169031384 
H  0.5821224931462822  -4.575311302581082  0.36340354797367447 
H  -0.10452303963364695  -3.3353805730621935  1.0334426798966871 
C  0.17624003533134247  -2.730436670780378  -1.631807449233485 
H  0.5353176052783343  -2.1618233530195097  -2.317192876293288 
H  -0.6570098574222627  -2.3690617972949966  -1.3215334621840895 
H  0.029579759876131817  -3.608944606099671  -1.9914518585941983 
C  1.6805196769234452  2.6259348635942192  2.3690979964916656 
H  1.5773595953476978  1.6567155171260988  2.4749279827961006 
C  0.6128632537095755  3.166610415610072  3.188425433736713 
H  -0.23442239588399283  2.9555211806531583  2.7896128801753637 
H  0.7083118448773804  4.119468738449751  3.2505894348887456 
H  0.6570512394151045  2.7827761367060546  4.0663312192613095 
C  3.0903205698287683  2.8889339701194015  2.887144338366958 
H  3.730775308152741  2.4974003553046087  2.287353169185178 
H  3.189021036386313  2.4971654554767477  3.7579130842668476 
H  3.2402594215445157  3.8357304130675676  2.9429087209678966 
C  1.4489905193369625  4.138662830871442  0.25599090119228785 
H  1.4167693781149027  3.9659677508935234  -0.7090606090442308 
C  0.1429961618955491  4.883493778894316  0.574305411069806 
H  -0.6030708330955317  4.295732980465759  0.4354831404658888 
H  0.06086820641869983  5.6474792844836035  5.054937625181832e-05 
H  0.1547375921794658  5.171069230690316  1.4909015644215617 
C  2.669082948049814  5.022445386327545  0.4614032427618237 
H  3.4628101817652768  4.525441099644493  0.25192420448572345 
H  2.703860387996845  5.312079970887206  1.3749298308246436 
H  2.6086935275542036  5.788133646959507  -0.11505541613975787 
C  0.3264858382596927  1.5251273962099343  -2.463614103717133 
H  -0.23801585985679585  0.9094807434829375  -1.951686489879336 
C  -0.49932751654109553  2.7793531312500335  -2.614432566100917 
H  -0.6350288918612625  3.178347670303266  -1.751703063567292 
H  -1.3496215492678776  2.5577616295998995  -3.001193479722916 
H  -0.038335484839690404  3.397530153060062  -3.1858211435277894 
C  0.5893409053789407  0.838291835424388  -3.7271916383900954 
H  1.1165239745288156  0.053857976836617194  -3.561636985472037 
H  1.0679261387236099  1.4260187498746972  -4.317274605823626 
H  -0.24320625986074385  0.5845405328185267  -4.133305223581753 
C  2.849164588579531  2.2249834928712087  -2.1468954167869407 
H  3.384176843316717  2.5024292314634877  -1.3757878306369642 
C  2.722460900270229  3.4632370293184547  -3.071749164600463 
H  2.223449642757476  4.149339593667479  -2.6205230172317187 
H  2.2674505825239044  3.2160836098492425  -3.8790373461001555 
H  3.597920365403464  3.7947392171702017  -3.2848934452137826 
C  3.649723723982822  1.124586719735858  -2.869504275295943 
H  3.728958061059858  0.3585775722960578  -2.297277182380432 
H  4.527194999239275  1.4554279784565824  -3.0828625482265255 
H  3.1967977900359186  0.8750404039320527  -3.6776106201708574 
C  6.330830090114896  -1.3493714777193946  2.0874398993675913 
H  6.504525670687477  -1.4001566651302075  1.1235638803541452 
C  6.241574673006612  0.1375496897590056  2.4504227119379527 
H  5.4890237649052525  0.5304112074651317  2.001793724219284 
H  7.045552303003041  0.5837210066624194  2.1743974478080665 
H  6.133156314173368  0.23071375180833908  3.399466155100811 
C  7.510460179566836  -1.988740203584665  2.782087136780216 
H  7.5536852356007484  -2.918146184225997  2.546901555421593 
H  7.406497782495076  -1.9027478387760093  3.7318994514980672 
H  8.318840949029962  -1.5496729330269108  2.506779429864524 
C  4.615168340904516  -2.0944949623874036  3.7360418054967823 
H  5.234484308107167  -1.5280511557263783  4.2409053481132615 
C  3.2747347013595114  -1.4590531171672787  3.886935109918719 
H  3.2789365640592276  -0.6001265562220823  3.460255139387282 
H  3.0763209434811545  -1.3521070981493188  4.819143595023232 
H  2.6091373131315727  -2.0168577864089796  3.4787138736309537 
C  4.684584330915175  -3.4606176190149425  4.389345875115041 
H  5.562661610228538  -3.8295592403190755  4.2673275911444035 
H  4.034472520805021  -4.041136674868391  3.986296051710602 
H  4.502482180754327  -3.3755763149010742  5.3281121737758825 
C  6.37448357203657  -3.6329396985413362  -0.6527553134536361 
H  6.393693841744331  -2.652779147387882  -0.6579388933276187 
C  6.015543463770323  -4.0711453753029065  -2.062789328480164 
H  5.12989031538743  -3.770247773400161  -2.275824827701912 
H  6.05077208478457  -5.028481778741702  -2.120385725175589 
H  6.639722415880481  -3.689529443420902  -2.685108640631491 
C  7.78392448811627  -4.098569113435668  -0.290286470564229 
H  7.9942520685303675  -3.8154137968048776  0.6027528709996806 
H  8.415522356717009  -3.717368419916519  -0.9046153349679273 
H  7.826490933453506  -5.056343028844883  -0.3398383090069886 
C  5.208902543525183  -5.503067290982556  0.4820709040635755 
H  5.86622142681644  -5.916936641080546  -0.11346873082080505 
C  5.4968100607075145  -6.005853052738375  1.8625457641481729 
H  6.363875720139495  -5.701592625623639  2.1401112964567037 
H  5.477937276941379  -6.966104521218042  1.8647874857102231 
H  4.8315365870575295  -5.673202991913654  2.468913212393068 
C  3.832559314623655  -5.976138942845913  0.02231418746201086 
H  3.6604296871134734  -5.647439909257623  -0.8632264534186473 
H  3.1614041963374717  -5.6433010089834985  0.6229014142496436 
H  3.8092362753656346  -6.93589153785435  0.01790712332167732 
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