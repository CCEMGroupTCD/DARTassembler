%chk=WEYAPANA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
N  1.5285951058079448  1.4681610955545714  2.5967510208140213e-15 
C  1.607606062363097  2.28217976736018  -1.1132004854005877 
C  0.33240283042452234  2.5753464071133045  -1.866721068125007 
H  -0.3771079931783792  1.9366961201392974  -1.5669528707818552 
C  -0.10113889187163383  4.009045969357126  -1.5092514187135437 
H  -0.17576011844497508  4.094763814204707  -0.5353894451159467 
H  -0.9690631799905038  4.199277164588807  -1.9208960499843837 
H  0.5661004323462989  4.644430177414782  -1.8433167452114598 
C  0.523593197108468  2.4348925763828912  -3.3751782319584827 
H  1.2013102779989553  3.0749953844776776  -3.6796907002392834 
H  -0.32421661540287383  2.616818814456315  -3.8307152948962138 
H  0.8174250226462465  1.5237632772546428  -3.5831149256207784 
N  2.6977102432643374  2.856734131696723  -1.5472396739579053 
C  3.852723675026237  2.6466930906054404  -0.8521505620399092 
C  5.054899453063212  3.2009522434081674  -1.3371406240367136 
H  5.0552268582728095  3.703924275725278  -2.143106806647769 
C  6.213631943760622  3.015227692590396  -0.6495031496174244 
H  7.014761970415026  3.4161112480681277  -0.9664306891911688 
C  6.247269452825193  2.2379207073723144  0.52366612150892 
H  7.070343693761766  2.1035075386829076  0.9778394628551069 
C  5.101504042577986  1.6764386687213586  1.0119089247876567 
H  5.127775210484323  1.1512790371696464  1.8030833411620222 
C  3.8689690907959324  1.8824043498142633  0.33115282981120464 
C  2.634482630318595  1.2901213417486284  0.7162934239517222 
C  2.568821324853068  0.34360972309339366  1.8636373692075714 
C  2.881626066290969  0.789288110827022  3.177530248570716 
C  3.1454330297885775  2.1546189393508675  3.4792944523665894 
H  3.081093340771635  2.8113016275921208  2.7958185190803837 
C  3.490348939487112  2.524531785866431  4.745889174927675 
H  3.676449913384757  3.4374619738807883  4.930590661999169 
C  3.5756776062558577  1.5769733168317288  5.787813992939084 
H  3.8538250242144896  1.8461538558043407  6.654837194082758 
C  3.2563245758439066  0.26839512618412953  5.543199089615172 
H  3.270773108145378  -0.36085069894859423  6.256106146297468 
C  2.9045876480805557  -0.15749203081332164  4.2392143914103935 
C  2.605459457900507  -1.5007074490471697  3.960022063330834 
H  2.627159689882981  -2.1390039396962557  4.66477508386852 
C  2.2835791379291495  -1.9088288853216189  2.696414860570254 
H  2.0842634186411892  -2.824537130193872  2.538141254691994 
C  2.242911983349118  -0.9907275890256236  1.6208697465568833 
P  1.528595105807944  -1.468161095554571  1.1102230246251565e-16 
C  2.7827496454742766  -1.1125987757378784  -1.2532852896930353 
C  4.141611427248517  -1.0665615661941472  -0.935987191099834 
H  4.426503367438482  -1.2238692354269718  -0.043587432238723425 
C  5.083430082348617  -0.7882919939809715  -1.9291087124081527 
C  6.554140515783918  -0.6996970057376218  -1.6077259358895768 
H  7.058728923948342  -1.2691073841567695  -2.225944133045441 
H  6.7058204787023366  -1.0018975451357925  -0.6876596758645029 
H  6.853050342011952  0.2294807317263341  -1.7011020997948938 
C  4.637386373449829  -0.5841154660048159  -3.2238893946765748 
H  5.276805385358156  -0.42143379181852314  -3.9071440032187477 
C  3.2895014856130556  -0.6089945546215886  -3.5603061155613243 
C  2.847486057694536  -0.33252727781591984  -4.965759290283887 
H  2.873683028854841  0.6323040777411022  -5.132023369495313 
H  1.9317522799022464  -0.6603228218179985  -5.0907217628213735 
H  3.447348307142345  -0.7880202833347706  -5.5925487953099555 
C  2.356173584378019  -0.8772996174602178  -2.558076594926744 
H  1.4291113481429245  -0.89947006598835  -2.7648722944924558 
C  1.4027767663384603  -3.277126639623016  0.0697526044409611 
C  0.46798927959237346  -3.869274461570045  0.9208223851207618 
H  -0.1788062595115445  -3.3306887693089937  1.363549640535409 
C  0.47967931183951995  -5.249398414680667  1.1243501650338885 
C  -0.4814368603910191  -5.883827799044761  2.095274523432842 
H  -0.296361272324553  -5.555858547946932  3.00076383702866 
H  -0.37428590581582943  -6.85754992962608  2.073241245197322 
H  -1.4004259823325744  -5.651541256353626  1.8453172154026483 
C  1.40005641748584  -6.025074206287059  0.4252694402832147 
H  1.4017125341802086  -6.965804170581814  0.5528518903125401 
C  2.3149871498206345  -5.461020358580932  -0.4530366748709863 
C  3.2921305311890965  -6.332798063272483  -1.1986687931777271 
H  3.760304622748942  -5.79676187665583  -1.8719361447028318 
H  2.808855663837006  -7.062448755722794  -1.6405947027967653 
H  3.943123276223174  -6.70751581166248  -0.5694406668558433 
C  2.3107040256233953  -4.080319016515808  -0.6283039817672739 
H  2.9302806805402586  -3.6806810566375896  -1.2268107316927939 
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


