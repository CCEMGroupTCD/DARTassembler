%chk=GAZOPUQU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.1242715630576094  -1.678127960705024  1.1400321509826735e-15 
P  2.3994586081341445  1.427349751382299  -0.9050776761128904 
C  0.924845466667348  -2.9517434030590586  -1.3481079434065744 
H  0.08395009483721916  -3.4463482600398656  -1.1300905608755463 
C  0.6811400405996808  -2.299229118904017  -2.7075942707666205 
H  -0.040520282745867275  -1.6408868256644398  -2.628689691426298 
H  0.42724760843975984  -2.985872111561086  -3.3579812659678496 
H  1.5011080199309743  -1.8524391406922223  -3.006711786537389 
C  2.015219364563029  -4.014146105205195  -1.4535105584613746 
H  2.8611760483993796  -3.5909721687542158  -1.7098295951738904 
H  1.7635091828920517  -4.675268728734983  -2.1309556123958453 
H  2.1215259072347004  -4.459170124521484  -0.5871217131864628 
C  0.9723086122904717  -2.746986827252024  1.5167198933877277 
H  1.4660116483231609  -3.596459053673364  1.3272679766624558 
C  -0.4833029315427939  -3.1264416244490594  1.784190097978017 
H  -1.0088909206111165  -2.3157217736050124  1.9457295952984357 
H  -0.5300369782610557  -3.708269013577238  2.5720675682344383 
H  -0.8448003156311683  -3.6006612034289986  1.0069526489256828 
C  1.5639745871447637  -2.1651920071539044  2.8012533579213756 
H  2.4322430255033494  -1.755389774485616  2.6044607283580294 
H  1.683182894536452  -2.8805515763345846  3.4603886750076462 
H  0.9567125813194705  -1.48536467770478  3.1628102115628445 
C  2.921042069473133  -1.2421132483986996  -0.010380651608003125 
C  3.823515578000868  -2.233689158558086  0.4062076233077487 
H  3.4882267374056024  -3.077410128087225  0.684539312033804 
C  5.195831903722437  -2.0134888468698438  0.4215664933987334 
H  5.784966328451527  -2.7040648649901353  0.7024944072774411 
C  5.700062750365033  -0.794809820068572  0.030126609216226993 
H  6.638518159021603  -0.643595677916152  0.03170748757102201 
C  4.83057762758207  0.21564838973779027  -0.36752268035191876 
H  5.182637860881153  1.0591954273577828  -0.6275311560828468 
C  3.4426181041107142  0.013342598433669868  -0.39132053760864305 
N  1.124271563057609  1.678127960705024  1.1267562259864795e-15 
C  0.9148010146860764  2.774357225145105  0.8847810741364807 
C  1.0384762224893143  2.5638517849928615  2.265584261647728 
C  0.7703943208095094  3.6026476371682974  3.1519313559030686 
H  0.8568315933379841  3.4468795008158235  4.086044467691519 
C  0.37931151765015236  4.8617479736601  2.7103063523998028 
C  0.28935812516355797  5.060723944395798  1.3419132618552658 
H  0.0454843077085183  5.921789809952873  1.0215830005118274 
C  0.5436216835217369  4.0455491842515965  0.41234074222206213 
C  1.4465052733790982  1.2084209696667227  2.7757552609159286 
H  0.7946489380876461  0.5389133181000433  2.4836383619968787 
H  1.4831214768503396  1.224830403657938  3.755164789213182 
H  2.330925954953453  0.9793870556864529  2.4211262515970784 
C  0.04331188781237283  5.969929875485223  3.6841284620359223 
H  0.1504016650865878  6.8369475905042405  3.2410515483442173 
H  0.645622906060962  5.920116147356361  4.455141230922581 
H  -0.8836177901465827  5.868603600623223  3.985780675920575 
C  0.3795455870272757  4.344852522114451  -1.052653249746282 
H  -0.4057114156182109  3.871018703347772  -1.3969455557730743 
H  1.176936881392949  4.0495929640818  -1.5379861221753184 
H  0.2595429545126021  5.309636872838102  -1.1771020536308234 
C  3.503100220788565  2.86122552626009  -0.9308455444035124 
C  3.9911819991140973  3.3485815647450243  0.29313223702020225 
H  3.7630721374485825  2.9096163086685327  1.1040156780699513 
C  4.803818222983246  4.469231136193626  0.3159593189127152 
H  5.134757954455987  4.799560317269372  1.1429586334605542 
C  5.1368420809117445  5.1104754990775545  -0.8697656437662612 
H  5.696247431914989  5.878655143410765  -0.8472928312789766 
C  4.663688516241429  4.643208619935294  -2.0816032392783184 
H  4.894760091405989  5.0882901348670915  -2.8883503013845258 
C  3.844295049717137  3.5170414209774576  -2.1144220907022637 
H  3.5166610349821363  3.1943799537479167  -2.9468073500346947 
C  1.9554418615570448  1.0832763353375399  -2.6242208269675484 
C  2.8909759861051922  0.5028154391429716  -3.4785513893349305 
H  3.731902025882607  0.2171083439590214  -3.138460076029454 
C  2.5951809701122297  0.3418358468081247  -4.8254775600530735 
H  3.236367172209242  -0.046619303969609  -5.410549094204814 
C  1.368214383690194  0.7477191292211307  -5.314290056682957 
H  1.1729916029688157  0.6500520181081356  -6.23837577851522 
C  0.42073530347937427  1.297240230468128  -4.463942850038958 
H  -0.4266652932759103  1.5569435815861803  -4.803821075344562 
C  0.7100223985941488  1.4689455930788229  -3.1146032133830426 
H  0.06256527405030932  1.8470507375359995  -2.530753522441552 
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

