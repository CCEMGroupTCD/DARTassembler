%chk=OPIKEPOL_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4656093297669741  1.389996148375959  1.7771262570476898e-16 
N  1.4656093297669743  -1.389996148375959  -2.220446049250313e-16 
C  2.656648922564785  0.4269133106543679  -1.0500027981268447 
H  3.4734127998896467  0.2629646765043674  -0.49781277281162284 
C  3.114194218675487  1.0934319704955873  -2.352508716766171 
H  2.326254689657808  1.3126434340352588  -2.90965319790568 
H  3.590960050857076  1.9357907794520193  -2.146133771005551 
C  4.040886577671766  0.14984143858481147  -3.1157911807769167 
H  4.286572083139951  0.5637203659992196  -3.9820628794056647 
H  4.871621441064173  0.01593044881877792  -2.5966381266933523 
C  3.3841221262037737  -1.1986031432148798  -3.3715145569141445 
H  4.019252329811634  -1.7906807406683387  -3.8469147073756336 
H  2.589155643485559  -1.0750368259983336  -3.9474600461786915 
C  2.96330051481234  -1.8557000215531052  -2.0533660357435926 
H  2.5221377693814997  -2.723322199074284  -2.2361822798723425 
H  3.762266520308353  -2.0267586646132365  -1.494289535883161 
C  2.0021806162792646  -0.9361215561786265  -1.3100506742291718 
H  1.2195798365693793  -0.7776107284634366  -1.9122158235815272 
C  2.1075564095997006  -2.275317184292352  0.6646988040144057 
H  2.7769274552491012  -2.747186466639886  0.1823032467970248 
C  1.9481969512257957  -2.6551397547482836  2.0665596344819095 
C  2.5880690369951775  -3.8290338312940726  2.4733808616860307 
H  3.05720652204539  -4.355350822067579  1.8374212518491368 
C  2.5437170825929623  -4.22817204160366  3.7943547952712686 
H  2.9595246661695223  -5.040998013330317  4.056821484827911 
C  1.900789600786059  -3.4557071781185957  4.734207415535961 
H  1.8817106639414058  -3.7296729840274336  5.643644626871455 
C  1.278483887937165  -2.2744098328572613  4.351809622872154 
H  0.8448017204255771  -1.7364237126268476  5.0041051402734915 
C  1.2879518373279581  -1.8739330073276426  3.0246709262213054 
H  0.8482041126766566  -1.0736992647380592  2.7653546986437068 
C  2.2571553581140624  1.9591952867508415  1.5200848614078424 
C  1.738819147895927  3.0593389646931715  2.2225371665932068 
H  1.026310629914871  3.564971578007366  1.851625624365795 
C  2.257905383932607  3.4054160747355584  3.445081674356424 
H  1.8911009386162614  4.143831800706611  3.918250327729435 
C  3.3090354585427466  2.6942042370948296  3.9996129974120604 
H  3.687504088369383  2.9636978111868655  4.829350744992743 
C  3.8006786982814296  1.5889141276257948  3.3345747212736114 
H  4.501552575193228  1.0790831177859572  3.72577785273901 
C  3.2830116409402486  1.2113287396279837  2.0992942662698493 
H  3.626632619595702  0.4467920743705731  1.6511040717040173 
C  0.9318976079466297  2.8094631571912525  -0.9699752310634355 
C  -0.07072082943248548  2.602675569804171  -1.9084577576616193 
H  -0.4802303839489195  1.7487696712973446  -1.9814772841939985 
C  -0.47805417065965283  3.6370377093538724  -2.7418569857575594 
H  -1.1473930214166272  3.4834521830562415  -3.399528992900369 
C  0.08696249156284575  4.883673437103614  -2.61479265673945 
H  -0.20375579375874686  5.593802793327571  -3.175001034049745 
C  1.0788034821454402  5.1104492021791454  -1.6744208463502905 
H  1.4593187501573435  5.975939744297261  -1.5814679963997516 
C  1.5170833252473515  4.070565435818487  -0.8663055194382432 
H  2.217379666224239  4.218610774102649  -0.24179900816175695 
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


