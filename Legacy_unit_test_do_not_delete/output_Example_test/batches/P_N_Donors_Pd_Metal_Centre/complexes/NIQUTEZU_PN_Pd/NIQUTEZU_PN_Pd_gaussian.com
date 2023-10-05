%chk=NIQUTEZU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  2.549113716696038  -0.9722539426212745  -1.1065135157547108 
P  1.5078446421962715  1.489464512836744  1.4037601005161293e-15 
N  1.5078446421962708  -1.4894645128367445  -2.220446049250313e-16 
C  1.4305939773011673  -2.8803733811885985  0.32256696159624726 
C  1.7463082303925521  -3.2982421250464924  1.620037174472676 
C  1.7522086226205325  -4.658002809201838  1.8985239554820146 
H  1.9720354575710837  -4.951262394264459  2.7745060776766026 
C  1.4457839364437088  -5.59031401100098  0.9324481123406545 
H  1.4672115417821572  -6.517814433531304  1.1414056054715236 
C  1.1123956066819973  -5.182650216552087  -0.3275575036286806 
H  0.8857050849556553  -5.82305395009478  -0.9907244776435454 
C  1.103204040122469  -3.831864909040351  -0.638345225578258 
H  0.8708395841060668  -3.5537093192767903  -1.5146980012705038 
C  2.0724359940834978  -2.3380244712024885  2.7106265055209393 
C  1.0616769048634116  -1.7852946458828869  3.485759528789301 
H  0.1542978515527409  -1.9886727407547606  3.2864534974220794 
C  1.3569138535602492  -0.9394421101950555  4.551224423562799 
H  0.6519454077150073  -0.5751344577981868  5.075180542686954 
C  2.6633981371766016  -0.6261412453686445  4.850561499052192 
H  2.8638998415785872  -0.053774050546164354  5.581767838482405 
C  3.6732801197272997  -1.1516426005565488  4.081385916356577 
H  4.577073663683006  -0.9232335165750439  4.271669547441223 
C  3.3875720959606563  -2.0107188719092783  3.025960563641336 
H  4.099520809377233  -2.37837280745492  2.5140812095525886 
C  3.051236877661303  0.6249766993287413  -0.4882650077595816 
H  3.5280507060767916  1.1369463289694464  -1.1907215625423966 
H  3.6540944318941606  0.5231935425677467  0.2916622057059879 
C  3.963330152839222  -2.077218223955845  -1.2559643394971214 
C  4.971747302319631  -2.0445354815471073  -0.3072132335781306 
H  4.96031061786428  -1.3791007253256706  0.3695519601642676 
C  5.999752268002626  -2.9724862102070606  -0.3397861626847578 
H  6.692580087135515  -2.9421260810396497  0.30890080002722975 
C  6.013081302279255  -3.944123105565995  -1.318801019784765 
H  6.71136974669983  -4.587859765108076  -1.3360420043936667 
C  5.017461933414527  -3.9816985425889335  -2.2706129018675476 
H  5.035407687820142  -4.650193203870693  -2.944818996275357 
C  3.993299856239356  -3.048935802673774  -2.2459213911355755 
H  3.3126622720905896  -3.0752055557714693  -2.9076328305986747 
C  1.8576960983626476  -0.7156034875715273  -2.7394272033025286 
C  2.6718277444701437  -0.2116121178340198  -3.758377134376892 
H  3.587392572811922  -0.0266495640936093  -3.589562920469027 
C  2.138385806858028  0.01505294003769131  -5.011933615249092 
H  2.6888717617818942  0.34277709915198296  -5.713339281724916 
C  0.7864791777510068  -0.23974266719514548  -5.24453715071524 
H  0.4203960470474728  -0.0888123240667611  -6.108820433520999 
C  -0.02567354681599232  -0.7144801878986707  -4.229528816555022 
H  -0.9489888260063086  -0.870480232208198  -4.392211671795316 
C  0.5082412299596384  -0.9589295264777827  -2.978682646153212 
H  -0.04509124202485992  -1.2935612217645238  -2.2841029070612384 
C  1.0846261333713647  2.4957949470847534  -1.445870226159419 
C  2.0603415062036614  3.2828543599789306  -2.032378343574204 
H  2.9425254017187132  3.284905389039155  -1.6817709024471046 
C  1.7470323581747715  4.068019968102543  -3.130536818285999 
H  2.412460102119688  4.613719662855249  -3.5282573995175994 
C  0.47098634700458875  4.052808951119768  -3.644956678263262 
H  0.25141877333404494  4.610007318618549  -4.383507610731116 
C  -0.4819438997678689  3.2365837241859756  -3.091770675454912 
H  -1.350676875997327  3.203301883918707  -3.4725747919242194 
C  -0.18363165324887687  2.4584363631638624  -1.9791731165365298 
H  -0.8511440591546791  1.9056505109975101  -1.5903125893464078 
C  2.027502854812155  2.6344275786247504  1.2927855429429265 
C  2.6126710876161097  2.091097195104092  2.435961038507993 
H  2.6864729879994114  1.147001358727406  2.5314378484464157 
C  3.084712778191192  2.930710300500092  3.4304687436598043 
H  3.5001260246497647  2.5635518752493525  4.202012141942872 
C  2.9534380679640657  4.299497151364437  3.3062352242068735 
H  3.2837229208678194  4.872648668330837  3.9889371390580637 
C  2.3471299094432743  4.8274411461510995  2.200885679818388 
H  2.2451331408718644  5.76862797062273  2.126518276730746 
C  1.8800269475889049  4.007674713666895  1.1886148518286357 
H  1.4573104624745241  4.385402195692287  0.42451844918345744 
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

-Pd 0
lanl2dz