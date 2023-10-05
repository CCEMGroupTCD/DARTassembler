%chk=DEYOBEHA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.374021636474478  -1.4805960092138564  9.655790115279554e-17 
N  1.374021636474477  1.4805960092138564  2.627684934986532e-16 
C  2.20609067053302  -1.5562740600141516  1.6006413730714453 
C  3.583420718189207  -1.5800114602106166  1.7213360525537542 
H  4.111979426067977  -1.5376990572255083  0.9573667188242132 
C  4.18245889719287  -1.6661217017664443  2.9629844312806877 
H  5.10932807386094  -1.6780695017671972  3.034278816878614 
C  3.4077409705695927  -1.7304841709389767  4.091121327552829 
H  3.8062889711443915  -1.799185237514192  4.928056883428685 
C  2.037128984688123  -1.694506564847573  3.982547085589055 
H  1.514602621841179  -1.7438330384740808  4.751001241104845 
C  1.4270026255179264  -1.5877351910114805  2.750826943034378 
H  0.5002160777066923  -1.537093848245866  2.690112865597168 
C  0.9169346760463541  -3.169655857030036  -0.3907978775553549 
C  0.30831111803377187  -3.458237406177169  -1.6121070702350715 
H  0.12486478988134397  -2.763100127949277  -2.201592361913595 
C  -0.026108333061692024  -4.7387425326273345  -1.9640162372993357 
H  -0.4097098931447183  -4.913228059666658  -2.7924609391085102 
C  0.2140528204206542  -5.772384323146944  -1.06830427508745 
H  -0.013058435786079903  -6.646250662748239  -1.2941029416380552 
C  0.7783674987385426  -5.512060083718085  0.13262900957111096 
H  0.9300546079798873  -6.211866187849183  0.7256205925728193 
C  1.1341032829489262  -4.226935064361663  0.4962569876701363 
H  1.5162193399094184  -4.064083821686746  1.3283645719236472 
C  2.646240499664511  -1.0585839234648469  -1.2195112355469686 
C  3.063562606485437  0.29082683582852004  -1.3264538456789434 
C  4.103065657407074  0.6033680166234774  -2.21186364200603 
H  4.408483160584693  1.4797151072042045  -2.271275893128488 
C  4.677245873989557  -0.36683066084203503  -2.9967566443505826 
H  5.35659528966692  -0.13929903985727177  -3.5891853597846377 
C  4.258908710220549  -1.6708952523879634  -2.911102639856678 
H  4.637303031368539  -2.321108861643306  -3.4560761989052287 
C  3.2599751518176783  -2.01129753083771  -1.9989646652482145 
H  3.0055609518618165  -2.902224621543834  -1.915514596903479 
C  2.514429252773003  1.3983384256177303  -0.5485486018843646 
H  3.0699301107870696  2.1368067236516697  -0.4446144903304755 
C  1.088855406433456  2.7239549494945727  0.7710938573664178 
H  0.25663779556163924  3.108863453538558  0.45372989963106347 
H  1.7952311876296325  3.368445073108145  0.6070331176872819 
C  0.9816226421012548  2.4777794683155254  2.264723512508546 
H  0.31634263132596496  1.7699924227753838  2.3937008183743234 
C  2.2517387095868604  2.0066902342101716  2.8972923456207784 
H  2.646541564334963  1.2165018285145484  2.473296827499465 
C  3.28450650171506  3.115251327030817  3.282356510438757 
C  2.3909151178227805  3.327452479840096  4.546220648461476 
H  2.881055224994434  3.55735759769918  5.362947972679317 
C  1.1908737553436892  4.178040639718145  4.274606509259196 
H  1.4683964924992154  5.103771525436093  4.183385253021567 
H  0.5837276686417288  4.123975874919475  5.028558563842082 
C  0.44627111079674986  3.744419327081428  2.992142962424083 
H  0.46866441329014463  4.483099141287266  2.365894727226344 
H  -0.4831199235444328  3.5889880996784913  3.2223620946313636 
C  1.997248184330293  1.858236031556534  4.424596376345875 
H  1.0749738662934814  1.6730798489827625  4.657642059301901 
H  2.6044338654415284  1.2466466416318678  4.870169329727327 
C  4.62229428005284  2.502163260231302  3.5787561434986106 
H  5.248550017837964  3.1932563697847303  3.804900068508805 
H  4.935811022833114  2.028566147557526  2.8047800670831395 
H  4.537921347682979  1.891123132222066  4.314733663702729 
C  3.479079145408422  4.3400084982748695  2.389191476521476 
H  4.1478112113686  4.911633678089707  2.7736765058247013 
H  2.6513824362642016  4.819600142843078  2.3186294319635676 
H  3.762336635247728  4.058457756272439  1.5163020080732563 
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

-Ni 0
lanl2dz


