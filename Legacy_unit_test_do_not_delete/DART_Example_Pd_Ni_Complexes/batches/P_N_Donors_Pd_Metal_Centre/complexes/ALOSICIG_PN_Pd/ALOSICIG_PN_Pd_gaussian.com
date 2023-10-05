%chk=ALOSICIG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.3892422349971942  -1.6006579936076284  3.963831245458204e-16 
N  1.3892422349971942  1.6006579936076293  -3.6255752253389427e-16 
C  3.458223440078572  0.2576482353800645  0.3781804626439879 
C  3.1218618238686062  -1.044894549103777  -0.024681855933597934 
C  4.111270892314016  -1.8986451725708082  -0.4872043615564426 
C  5.441595425293871  -1.494873459574616  -0.567787304764827 
C  5.7735800904135415  -0.21089259157677812  -0.1576705009870518 
C  4.8199824800742155  0.6341310020963885  0.2853433291479085 
C  2.487449396176765  1.216382750266185  0.9414181528358277 
C  2.022623514786467  0.8519159857669556  2.317896005024232 
C  1.9506770097402706  1.9762013060457706  -1.3247685511787464 
C  0.6771823818225737  2.7746770829476164  0.600151908571707 
C  1.1619473449294841  -2.549032437781438  1.5337069200427957 
C  -0.005441307042645782  -2.420055460298135  2.285215359467121 
C  -0.19796575299426955  -3.2151066799104826  3.397945670865173 
C  0.7320323163710917  -4.141437630008162  3.7315969586910236 
C  1.8881685657737965  -4.272451332396487  3.0302220927287578 
C  2.101859931230819  -3.4823203893400607  1.9153464614792957 
C  1.2809869980382964  -2.8649778190056843  -1.300224154819635 
C  1.5322124441248333  -2.3403398961304176  -2.7149245666099615 
C  1.2744196145465918  -3.3804217541216866  -3.7868764665082404 
C  -0.18531795696907838  -3.6700387315468084  -3.9312335193068515 
H  3.911658990472074  -2.802193207676254  -0.8039276901890754 
H  6.065888416781003  -2.079293194293605  -0.9171215923047589 
H  6.616229825893445  0.09234416914187735  -0.19369757753394415 
H  4.959562732885571  1.5013875869157725  0.5346636033185908 
H  -0.6641590793239083  -1.8177672551662152  2.0317470489433385 
H  -1.0438129854234652  -3.1178855148511935  3.849959713187306 
H  0.6421282590875319  -4.612292892662572  4.364759498734742 
H  2.587800218629434  -4.959589730725451  3.1246778594367437 
H  2.7990220741817105  -3.5929384168681855  1.3578847153281053 
H  1.9459115565298708  -3.6331153755044623  -1.0691799040561187 
H  0.5722646587088476  -3.161656411119978  -1.23959356257041 
H  0.917324446918862  -1.6735009965886787  -2.7767272603214104 
H  2.483256973828781  -2.018996569164791  -2.798821162100243 
H  1.6209837853876603  -3.092899546480862  -4.567320869142176 
H  1.902055318308932  -4.175142747038445  -3.627062772538795 
H  2.8669924313485495  1.8847350649220052  1.0418734062913086 
H  -0.46955124106143864  -4.323644295120094  -3.185439326864908 
H  -0.2706199303534369  -4.021966121563696  -4.622818001125708 
H  -0.7326611204135463  -2.9141858449016333  -3.814907704541751 
H  1.547694049747678  1.5306302703540406  2.733267288498926 
H  2.589037126143081  0.592952408503486  2.8438333078663938 
H  1.5712211500754893  0.1292311151577002  2.258474519635333 
H  1.1285188473666934  2.1961669710448537  -1.9358417892312985 
H  2.6173084049464928  1.1729318317828135  -1.955806834203404 
H  2.669032032503658  2.708546000965872  -1.1456509666343508 
H  -0.017742716235161993  3.0487779381909936  -0.006962981086192392 
H  1.2869190376685018  3.510360411118131  0.623542957496098 
H  0.2042382611265663  2.473286531732678  1.4062064131416192 
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


