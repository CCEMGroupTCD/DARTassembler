%chk=EFOZEPIT_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4628804069369437  1.5336495411273074  -9.914251097415421e-17 
N  1.4628804069369437  -1.5336495411273074  1.1102230246251565e-16 
C  1.3095398017659003  -2.4069167798835363  -1.1588030240385456 
C  1.5533130165156588  -2.1080727234616816  -2.467519608984235 
H  1.8921619941728132  -1.2978693903803489  -2.748195951778919 
C  1.3000724412709528  -3.0986081354068142  -3.3953539954629983 
H  1.4628611571415076  -2.956812556189834  -4.227549031316822 
C  0.8154421820763914  -4.336023610134536  -3.021324645307985 
H  0.6049149885473781  -4.9214021343181615  -3.6428889494707937 
C  0.5744272401419621  -4.628380486967313  -1.6974892498252463 
H  0.3168880770459208  -5.535921921744466  -1.4996987201368452 
C  0.8351545800208928  -3.6563588842325476  -0.7323899036143934 
C  0.7207782196116072  -3.639538379085695  0.7038840357113116 
C  0.257107485530641  -4.562551829707199  1.634570637274166 
H  -0.0441810841614807  -5.333322115421303  1.3123661456340332 
C  0.23094969150250577  -4.194684214016263  2.9545830149197054 
H  -0.08522195235745511  -4.80150381330746  3.5846507934443648 
C  0.6573209245652154  -2.954112050051415  3.3812837386268173 
H  0.6284836474976165  -2.7329984573392654  4.283681943727707 
C  1.1275949641741851  -2.041689304518939  2.453297752145963 
H  1.3676951122610528  -1.2460662744340554  2.710300034017554 
C  1.134272668555419  -2.38153993924536  1.141088251212421 
C  2.7887888852690703  -0.9019489380756617  0.06607976385556513 
C  3.900319494865915  -1.7251906345273507  0.14872529097359985 
H  3.733462796460814  -2.6584947086431887  0.2300503021087341 
C  5.166843977591174  -1.196031052603912  0.13557456794596634 
H  5.849942878666109  -1.7292917018385796  0.3350191249739093 
C  5.346634284578047  0.16930982146621143  0.02046087151065673 
H  6.203235735899686  0.5303359947554571  -0.012846121655930645 
C  4.238016770288142  0.9940346175042163  -0.04580282752411198 
H  4.396516970685016  1.9311582213384906  -0.13882442507445586 
C  2.9455554165659628  0.48025923801248593  -0.0021225280429941484 
C  1.5788330722284831  2.6147779968377516  1.475079432323899 
H  0.7663098435921216  3.1632183298109147  1.4738571006816668 
C  1.5205566081487907  1.7749304107632748  2.7580053911770177 
H  0.7228861308462599  1.200127676266304  2.7122882603095015 
H  2.3547823042693903  1.2560909245172138  2.752731356286031 
C  1.5084534564875896  2.6585041501743922  3.991105853727818 
H  0.6977624323574911  3.0497706546190257  4.023838649126023 
H  1.523061004552194  2.0777323922763706  4.787413146654114 
C  2.6536625634227495  3.6353981339842405  4.023568998814278 
H  3.52438911173462  3.0952693185111104  4.0355594433643995 
H  2.5012782103189277  4.133504392626929  4.743896357749158 
C  2.677853568239599  4.473416772228864  2.7796329927268695 
H  3.361114100947188  5.131014879198391  2.727740949461498 
H  1.9120104771138835  4.959143018690908  2.662574232811575 
C  2.7536991673206312  3.600489581825284  1.5330073170109761 
H  3.589709163438589  3.1062602193316535  1.5357337872546144 
H  2.7415626702981424  4.164352799470402  0.7432761303234161 
C  1.7343425286817682  2.622461264549187  -1.456573169256147 
H  2.692114542585416  2.9970096568616653  -1.3821122318518722 
C  0.7736911023856639  3.810247027299782  -1.4986585463157835 
H  0.8978191170116042  4.359340697279304  -0.7096664629756184 
H  -0.14172806281329486  3.4901201457443523  -1.5054225460232957 
C  1.0386840818505552  4.641530270838517  -2.7540480873929947 
H  0.40671243197759366  5.376615855757306  -2.792336944799008 
H  1.9331841964179297  5.0167664897839455  -2.709129215740869 
C  0.9098393677139935  3.807305628447551  -3.9962784869834724 
H  1.0780896816523988  4.360148861052554  -4.776289161283846 
H  0.0050547067677593205  3.462975052368328  -4.062268650594277 
C  1.8970853053050356  2.6421333892366987  -3.9755210381425576 
H  1.779438685019861  2.1003578937807714  -4.770690972973443 
H  2.805906189969714  2.9834396009643434  -3.9749132277959505 
C  1.669584020772128  1.791350159424686  -2.737141592493288 
H  2.3417979715433326  1.0940387846148236  -2.7006066646931424 
H  0.8001540889454035  1.3656746465430851  -2.7979063962596817 
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

-Pd 0
lanl2dz


