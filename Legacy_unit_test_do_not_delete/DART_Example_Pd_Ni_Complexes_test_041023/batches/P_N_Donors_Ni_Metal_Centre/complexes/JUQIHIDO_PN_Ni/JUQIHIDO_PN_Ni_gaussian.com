%chk=JUQIHIDO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
C  2.191154644164575  1.2568671098551996  1.049826712774887 
C  1.8418777311715382  2.310911081558607  -1.0110189737736888 
C  3.875321201690044  2.610853386786566  -0.024206921891354205 
C  5.161403389346568  3.1694953927064926  -0.06197259167816138 
H  5.417089523698674  3.711863591654077  -0.773754178592604 
C  6.040679882250881  2.906396844843351  0.968776637348096 
H  6.899251865772371  3.263449019609644  0.9391978678222684 
C  5.6589796194677655  2.1009332474667435  2.0714762984945243 
H  6.267379525816955  1.9324705596763934  2.753664085527393 
C  4.4167201315139994  1.5764624340603846  2.1344891580846688 
H  4.1743199945368525  1.0621812240086907  2.8688398240346973 
C  3.47493036358516  1.7980143345727433  1.1023163632465174 
C  0.88238516085586  2.60399552585093  -2.137151046855221 
H  0.09457317527069065  2.0326513792019654  -2.0184202527943906 
C  0.42236366679779325  4.071061312545891  -2.048770661510484 
H  -0.19471275466538773  4.256290819524672  -2.760628835459912 
H  1.1826107330586264  4.650812885320514  -2.1263421260264335 
H  -0.009969767047493638  4.220591497426483  -1.2043973401714445 
C  1.4974347024130452  2.2854621814401312  -3.5177746885218353 
H  0.8589384134361617  2.48523212600605  -4.206711803366808 
H  1.7293243345777836  1.354179217504128  -3.557186610285103 
H  2.2857955836784556  2.817266906668774  -3.6474518457813945 
C  1.62846575696085  0.34388782858623723  2.1059797784765677 
C  1.2876250437733725  -0.9473997654124473  1.7653770071740666 
C  0.6714262590546  -1.7904800597094257  2.7306810085535376 
H  0.4435467110713269  -2.6637656677207557  2.508665314398331 
C  0.41131699634605723  -1.3106960337580678  4.001724560087171 
H  0.012133818029671728  -1.8709864984410334  4.626956041487911 
C  0.7424014081302117  0.015028744099210511  4.370579498577623 
C  0.469924197373462  0.5232076723118111  5.667654826590217 
H  0.0931634788757687  -0.03521278651392801  6.308492280366107 
C  0.7481651469755601  1.8118844994610346  5.987046669211705 
H  0.5653092625880938  2.1256614202138673  6.8433267983445285 
C  1.3150252675103824  2.6782519209540823  5.026371482313056 
H  1.4806353697572707  3.567516912480855  5.242308901207192 
C  1.613455644714205  2.2174625849817486  3.7998189853459508 
H  2.005579079051023  2.794117204141765  3.1840861598076806 
C  1.346699288054577  0.8635860485562514  3.414831166969581 
C  3.0972184592789267  -1.2277895685866962  -0.4476183924675693 
C  3.387743004244701  -0.5358109630702212  -1.6280849816940024 
H  2.706787222685576  -0.32466358208950896  -2.223631543321388 
C  4.687076951305439  -0.1681385711063219  -1.902983119561822 
H  4.8726775541374385  0.32878444055053135  -2.6665703357423127 
C  5.711860391118087  -0.5278767234896315  -1.0549400350397988 
H  6.584739140792298  -0.27720481633559646  -1.2584407557534774 
C  5.467220826158416  -1.256105708162668  0.08910487628949981 
H  6.1638702320862535  -1.525453501395137  0.6423647713670112 
C  4.139220489721188  -1.5777587398030615  0.3980441157189535 
H  3.9527082385654895  -2.0355034953384266  1.1854131136790707 
C  1.0562405220086102  -3.274328188043418  -0.06391753418331557 
C  2.077376146975687  -4.172011419088402  -0.40304729124358263 
H  2.937065705641234  -3.85895155304413  -0.5706669887549134 
C  1.7995944575442533  -5.542595995327493  -0.49000695081131257 
H  2.475750537549092  -6.135560520884179  -0.7246873226020494 
C  0.5472584163982588  -6.014308276532481  -0.23041855775954861 
H  0.3770050727274793  -6.9271319274578245  -0.2710514419820991 
C  -0.4645798534202028  -5.139686788579117  0.09000464449186674 
H  -1.3176213963141061  -5.470485028958809  0.2554121118531808 
C  -0.23430609188952212  -3.7785784537731066  0.17112841210729124 
H  -0.9311887749683587  -3.2000174744252554  0.38082707495684165 
N  1.3630048954424192  1.490743993782971  -1.1895134953576015e-15 
N  3.0242998865843553  2.864094732433041  -1.056691227869393 
P  1.3630048954424188  -1.4907439937829705  -1.1102230246251565e-16 
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

-Ni 0
lanl2dz


