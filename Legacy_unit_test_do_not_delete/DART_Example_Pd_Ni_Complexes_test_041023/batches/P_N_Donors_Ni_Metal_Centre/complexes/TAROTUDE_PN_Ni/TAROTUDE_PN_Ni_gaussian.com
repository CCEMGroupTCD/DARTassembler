%chk=TAROTUDE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3145413972180555  -1.533649541127308  1.8781790015553647e-16 
N  1.3145413972180555  1.533649541127308  -1.8781790015553647e-16 
C  1.1612007920470164  2.4069167798835363  1.158803024038544 
C  1.4049740067967764  2.108072723461682  2.46751960898425 
H  1.7438229844539337  1.2978693903803495  2.7481959517789414 
C  1.1517334315520706  3.0986081354068156  3.395353995463014 
H  1.3145221474226263  2.9568125561898353  4.227549031316812 
C  0.6671031723575067  4.336023610134537  3.021324645308001 
H  0.4565759788284838  4.921402134318163  3.64288894947081 
C  0.426088230423076  4.6283804869673135  1.6974892498252465 
H  0.16854906732703512  5.535921921744467  1.4996987201368432 
C  0.6868155703020056  3.6563588842325485  0.7323899036143936 
C  0.5724392098927149  3.639538379085695  -0.703884035711314 
C  0.10876847581174665  4.5625518297072  -1.6345706372741717 
H  -0.1925200938803706  5.3333221154213035  -1.312366145634035 
C  0.08261068178360764  4.194684214016264  -2.954583014919731 
H  -0.23356096207634702  4.801503813307459  -3.5846507934443816 
C  0.5089819148463168  2.9541120500514153  -3.3812837386268306 
H  0.48014463777872707  2.732998457339264  -4.283681943727702 
C  0.9792559544552932  2.0416893045189375  -2.4532977521459736 
H  1.2193561025421626  1.2460662744340545  -2.7103000340175707 
C  0.9859336588365284  2.381539939245359  -1.1410882512124185 
C  2.6404498755501806  0.9019489380756613  -0.06607976385556596 
C  3.7519804851470346  1.725190634527351  -0.14872529097360415 
H  3.5851237867419385  2.6584947086431896  -0.23005030210873786 
C  5.0185049678723  1.1960310526039122  -0.13557456794597844 
H  5.701603868947204  1.72929170183858  -0.3350191249739141 
C  5.198295274859168  -0.1693098214662112  -0.02046087151066215 
H  6.054896726180775  -0.5303359947554567  0.012846121655925241 
C  4.08967776056927  -0.9940346175042167  0.04580282752410339 
H  4.248177960966149  -1.931158221338491  0.1388244250744506 
C  2.7972164068470766  -0.48025923801248593  0.0021225280429931076 
C  1.4304940625095939  -2.614777996837753  -1.4750794323239014 
H  0.6179708338732276  -3.163218329810917  -1.4738571006816625 
C  1.3722175984299068  -1.7749304107632762  -2.758005391177042 
H  0.5745471211273665  -1.2001276762663056  -2.7122882603095158 
H  2.206443294550509  -1.2560909245172147  -2.7527313562860525 
C  1.3601144467686919  -2.658504150174394  -3.9911058537278303 
H  0.5494234226386003  -3.0497706546190284  -4.023838649126008 
H  1.3747219948333087  -2.077732392276373  -4.787413146654094 
C  2.505323553703854  -3.6353981339842423  -4.02356899881427 
H  3.3760501020157117  -3.095269318511112  -4.035559443364385 
H  2.3529392006000327  -4.13350439262693  -4.743896357749136 
C  2.5295145585207126  -4.473416772228865  -2.77963299272689 
H  3.212775091228302  -5.131014879198392  -2.7277409494615106 
H  1.7636714673949927  -4.9591430186909085  -2.662574232811592 
C  2.60536015760174  -3.6004895818252844  -1.5330073170109741 
H  3.441370153719712  -3.106260219331655  -1.5357337872546237 
H  2.593223660579255  -4.164352799470402  -0.7432761303234187 
C  1.5860035189628838  -2.6224612645491883  1.4565731692561505 
H  2.543775532866529  -2.9970096568616658  1.382112231851873 
C  0.625352092666778  -3.8102470272997815  1.4986585463157827 
H  0.7494801072927167  -4.359340697279304  0.7096664629756206 
H  -0.29006707253219277  -3.490120145744352  1.5054225460233077 
C  0.8903450721316646  -4.641530270838516  2.754048087393008 
H  0.2583734222587031  -5.376615855757306  2.792336944799032 
H  1.784845186699045  -5.0167664897839455  2.7091292157408944 
C  0.7615003579951096  -3.807305628447551  3.9962784869834715 
H  0.9297506719335221  -4.3601488610525525  4.77628916128383 
H  -0.14328430295112393  -3.462975052368329  4.062268650594285 
C  1.7487462955861546  -2.6421333892366987  3.9755210381425696 
H  1.6310996753009817  -2.1003578937807705  4.770690972973424 
H  2.657567180250833  -2.9834396009643425  3.9749132277959465 
C  1.521245011053242  -1.791350159424686  2.737141592493303 
H  2.1934589618244518  -1.094038784614823  2.700606664693154 
H  0.6518150792265158  -1.3656746465430851  2.7979063962597035 
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


