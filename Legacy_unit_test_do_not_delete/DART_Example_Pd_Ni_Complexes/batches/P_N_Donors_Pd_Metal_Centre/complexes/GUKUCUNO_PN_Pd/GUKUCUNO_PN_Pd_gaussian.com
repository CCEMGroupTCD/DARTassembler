%chk=GUKUCUNO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.4646129147320808  1.5319951076945384  -1.1090949837783824e-15 
N  1.4646129147320806  -1.5319951076945384  -3.0531133177191805e-16 
C  1.0221061858462708  3.2033173117537728  -0.5717011939023521 
C  2.328821383233298  1.8273345168920612  1.5697031236337995 
C  2.809443996776337  1.0653188197552748  -1.1437933349355445 
C  3.206083273553464  -0.2802058948125512  -1.219179849677219 
C  1.3484533936414511  4.379232705315982  0.10589394860365958 
H  1.8136270782491186  4.3354758619270495  0.910554679474242 
C  1.2312626139765284  -2.794884297276575  0.6601984716538879 
C  4.271078046029421  -0.6230782364811084  -2.0692294140660437 
H  4.56937193521124  -1.5035934504167412  -2.091105284441085 
C  3.463453408521655  2.0037233047170164  -1.9329539281941448 
H  3.2155435627172717  2.898491364428939  -1.878557993855191 
C  -0.008731031177025894  4.539728870587695  -2.309172277809889 
H  -0.4542398108346095  4.591132853279593  -3.1233357858421837 
C  4.883765300984901  0.3144490002484095  -2.8695062024370457 
H  5.562513205151889  0.061963774899294055  -3.451982943233894 
C  1.5425292210366948  2.1149103308011497  2.6875417823510466 
H  0.6197784772006995  2.190559247313949  2.592641328292164 
C  1.1992952857731056  -2.8284260101331706  2.057980104608846 
C  2.1159189237125933  2.2885982744920605  3.9319022152970917 
H  1.5841238942423426  2.4969620795699856  4.665408664733065 
C  0.33796662652234  3.3051173104722387  -1.7750168605621923 
H  0.1056973944081947  2.530155827868915  -2.234246905631499 
C  1.0321251026874363  -3.946707864556424  -0.11013375632071958 
C  4.481947327204278  1.628122790648841  -2.8027515893994557 
H  4.891938312529307  2.2671917831558397  -3.3397279434117024 
C  4.279367593979751  1.872057029161116  2.9862274603958263 
H  5.200075031397992  1.7915324639190993  3.087200224423101 
C  0.31653982553589666  5.689154838827212  -1.6201420220809422 
H  0.08366979688830556  6.52031163796099  -1.9685228266775237 
C  0.9909667853445568  -4.048394801556134  2.6744229031418403 
H  0.9725951561748297  -4.096973507601742  3.602945731380646 
C  2.647109911373113  -1.4116619688521932  -0.48245629712060767 
H  3.218773212504768  -2.132836047314258  -0.3472636246355266 
C  3.70141783268937  1.71085622064702  1.729706134369566 
H  4.238607368417039  1.5238474968194735  0.9937117526538469 
C  0.9869776664285895  5.608767460707366  -0.411559231810373 
H  1.1941916249333548  6.385582539308313  0.054790324916120944 
C  3.4884835837642663  2.1515042400887365  4.084749785471778 
H  3.874797779725926  2.247702197418888  4.924788543508574 
C  1.4320634800209164  -1.5770300165225608  2.873204303784705 
H  1.3269564667585847  -1.7777334120201576  3.8056542181494706 
H  0.7950292855144283  -0.9056652681404297  2.6175919612437677 
H  2.3214922126072413  -1.2515851064860652  2.7132125534309717 
C  0.8080944873036581  -5.204132460969752  1.9238609235409654 
H  0.6793259049978415  -6.0196183876056875  2.3535546350026086 
C  0.8181442728832702  -5.151424590301995  0.5586045072661626 
H  0.6791522974699782  -5.928754311538574  0.06822442482706445 
C  1.0735757950164095  -3.9033516360238436  -1.6184417711349737 
H  0.5758290200650396  -3.143528961556075  -1.9303446107984852 
H  0.6860406514496664  -4.7070170223009695  -1.9729531805433995 
H  1.985000173293988  -3.831247136835047  -1.9122351575106433 
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


