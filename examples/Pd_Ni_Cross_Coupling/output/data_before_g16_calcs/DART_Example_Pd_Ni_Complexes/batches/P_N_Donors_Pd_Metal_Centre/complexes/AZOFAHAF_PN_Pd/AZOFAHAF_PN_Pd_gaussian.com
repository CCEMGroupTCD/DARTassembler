%chk=AZOFAHAF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5118602506514947  -1.4853883608336238  2.0324886605581075e-16 
N  1.5118602506514947  1.4853883608336238  -2.374187613898209e-16 
C  1.5594901162008972  2.4461630798431884  -1.2108649409491872 
C  0.8966686576077649  3.7580244161649574  -0.7908796456202334 
C  0.6223495098653467  3.566517711588144  0.6720210096821182 
C  0.07140957572292317  4.452866251616506  1.6329459362797991 
C  0.010001965573489091  4.02200125235071  2.9812632423407575 
C  0.49071802925992647  2.7726952735704873  3.396582790999157 
C  0.9560233519404274  1.8553658694035573  2.4065896794392967 
C  1.017912280369171  2.3121192565239737  1.122155545708811 
C  1.2457701500575662  -2.5698957007833796  1.4227750025780912 
C  1.8357336342348831  -2.2919994004728386  2.6681275685653674 
C  1.5780321373742847  -3.11740666925736  3.745858016847093 
C  0.745473831820739  -4.219553557685398  3.6135673362145235 
C  0.16350333239268866  -4.487256035693799  2.4154619731251015 
C  0.40203728318025767  -3.6671613102198855  1.3216728181360797 
C  1.9202012453973232  -2.554327078724822  -1.3988574549978836 
C  2.657272664024389  -3.716938517951946  -1.2120744653512576 
C  3.029900757225477  -4.496199905425313  -2.3465784204569387 
C  2.6473132887318314  -4.090412081921886  -3.6052957522070126 
C  1.9520543589090378  -2.946341982270455  -3.786291297967741 
C  1.5653466250567065  -2.158243252761747  -2.6995008116615886 
C  2.9649557044765107  -0.49095486859121945  0.345114649905925 
C  4.223478544632326  -1.1031232778614282  0.6178923240249311 
C  5.312652310344765  -0.3174260549849647  0.8291173747362601 
C  5.176290069031851  1.064871456655482  0.8169379471328587 
C  3.981771509731796  1.716979902191159  0.5758017890213203 
C  2.854494234313834  0.8952794720659266  0.30452516367682186 
C  3.9977644360279774  3.2236257760985216  0.5873869927918587 
H  1.0887952382392507  2.0666554205239547  -1.9437404818009176 
H  2.461317480153208  2.6053929528087068  -1.4625946237087486 
H  0.0905015708675061  3.8990594818244366  -1.271843041848221 
H  1.483450482262982  4.49256610827009  -0.9319266335839368 
H  -0.24196373059667176  5.315485179289331  1.381046876124933 
H  -0.3791142252773714  4.596102051498899  3.6266525221053847 
H  0.5087936577169967  2.5473904483649226  4.321713093926811 
H  1.2067803679889304  0.9674191513750239  2.6347393558467216 
H  2.407999060603224  -1.5395836218997734  2.769283576482641 
H  1.9738754445158584  -2.9301515665959257  4.587868184670786 
H  0.5794977484854704  -4.783303885263057  4.3623550483978875 
H  -0.4088783200135291  -5.241429920139281  2.326677931714788 
H  -0.015397419834299209  -3.86212056305628  0.4899581684799337 
H  2.912560793754672  -3.988709558202621  -0.3362873493807762 
H  3.5464369000854803  -5.2864284180152765  -2.2353681645936483 
H  2.8723668841129273  -4.624753100693948  -4.357131847951085 
H  1.7230180009082412  -2.6734400189263368  -4.668582681308115 
H  1.0643353166974645  -1.36359628626905  -2.8395058027873326 
H  4.299017232589503  -2.0511631354039292  0.6580189377931269 
H  6.164522341804954  -0.7107610207121143  0.9827804940475017 
H  5.9489063167055685  1.5941307606389856  0.9800027286692374 
H  4.597178163613654  3.5365126528606154  -0.07910723488901143 
H  4.285375430342819  3.5248699333833375  1.4408654598025965 
H  3.1252873288932532  3.547745159752374  0.41157998876272844 
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

