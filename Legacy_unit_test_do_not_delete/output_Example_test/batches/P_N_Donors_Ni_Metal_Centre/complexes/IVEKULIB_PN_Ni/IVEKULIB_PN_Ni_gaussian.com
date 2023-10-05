%chk=IVEKULIB_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3688394418265437  -1.4853883608336234  -2.6073782157449317e-15 
N  1.3688394418265428  1.4853883608336238  -1.5374156464789625e-17 
C  1.0603654618030038  2.446163079843188  -1.1718815268449985 
C  0.5492979653498476  3.758024416164957  -0.5764573517116978 
C  0.7146760181428051  3.5665177115881437  0.9027246593189184 
C  0.46875679777238  4.452866251616506  1.9827409844927744 
C  0.8042422374143571  4.0220012523507105  3.2900970844884716 
C  1.385380979831465  2.7726952735704864  3.5467213690873085 
C  1.5409086991952228  1.8553658694035569  2.4639441377233235 
C  1.2245611785102277  2.3121192565239737  1.2175390954074117 
C  1.5303553660545397  -2.569895700783378  1.4384037179512381 
C  2.458646104516364  -2.2919994004728377  2.4568516700625485 
C  2.527302825845509  -3.117406669257359  3.5628350493472687 
C  1.6924453068126537  -4.219553557685396  3.6797413333645284 
C  0.7856120429355234  -4.4872560356937985  2.704139183513415 
C  0.6939301946809403  -3.6671613102198846  1.588402834750615 
C  1.3503515383748144  -2.5543270787248202  -1.457121389870843 
C  2.1098265025884895  -3.716938517951945  -1.4939987557051804 
C  2.1344756649531957  -4.4961999054253114  -2.6878761944683998 
C  1.400592117213034  -4.090412081921886  -3.7797358148650484 
C  0.6827947197364386  -2.946341982270454  -3.7495486775481703 
C  0.6307310618435564  -2.1582432527617468  -2.5971833673823306 
C  2.8593432937054617  -0.49095486859121906  -0.09480921390212538 
C  4.142627144827054  -1.1031232778614273  -0.20190709503044393 
C  5.245965425204646  -0.3174260549849642  -0.3183551650773051 
C  5.112000645292367  1.0648714566554824  -0.2901339486668694 
C  3.8991754763302513  1.716979902191159  -0.1714901759980772 
C  2.7418412471555946  0.8952794720659268  -0.1013293242874721 
C  3.917856773585009  3.2236257760985216  -0.16508706965823028 
H  0.39596564017241764  2.066655420523955  -1.7351160282165838 
H  1.8491886225304808  2.6053929528087068  -1.6762806234493839 
H  -0.36226354186514875  3.899059481824436  -0.8007044895270428 
H  1.0692020666633408  4.49256610827009  -0.8828996594904074 
H  0.09542825691155254  5.315485179289331  1.8334702030273164 
H  0.6208221374160239  4.596102051498899  4.0210524862292525 
H  1.673148713049604  2.5473904483649226  4.4261430753566104 
H  1.8474133360566847  0.9674191513750241  2.6088105020898706 
H  3.0354814079033496  -1.5395836218997727  2.3862734232668776 
H  3.1520296112416366  -2.930151566595925  4.252319994894226 
H  1.7526459288295113  -4.783303885263056  4.444337294160087 
H  0.2122828049741694  -5.24142992013928  2.786582781925133 
H  0.05156557292126451  -3.862120563056279  0.9150762554111094 
H  2.6100151272609686  -3.98870955820262  -0.7311183971761814 
H  2.6609563670864116  -5.286428418015276  -2.732545850513348 
H  1.395996339837769  -4.6247531006939475  -4.564519552265337 
H  0.20580912548242214  -2.673440018926336  -4.526324373149585 
H  0.11067806694321791  -1.363596286269049  -2.584589275944091 
H  4.2265970378518585  -2.051163135403929  -0.18561919841378496 
H  6.105539535898437  -0.7107610207121136  -0.42046908656407095 
H  5.898532765635094  1.5941307606389858  -0.3600854519824091 
H  4.296214918529799  3.536512652860615  -0.9777102827797919 
H  4.442433489659625  3.524869933383337  0.5670091307513456 
H  3.0321017730930038  3.5477451597523735  -0.07812452459836054 
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