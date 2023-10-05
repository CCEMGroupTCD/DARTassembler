%chk=KOLIFUPI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5306830844430201  -1.465984070513728  -2.5287817477390755e-15 
N  1.5306830844430208  1.4659840705137284  -1.240201187243067e-16 
C  2.7048649880898403  0.9842149816091661  0.1876607367683344 
C  4.056436154412422  1.533119305459024  -0.16306057085979256 
C  4.920173168302533  1.555721689527487  1.1275151933131913 
C  5.224548172947795  0.05459638109979754  1.3892029938409716 
C  4.415090024591143  -0.6588551331676662  0.30490225011323513 
C  2.91912049676923  -0.4371752278149998  0.6505295513947538 
C  4.611017575325333  0.26953119269674886  -0.9031264725632191 
C  1.3149064657813498  2.8417070306501575  -0.3895034374590288 
C  1.3188947973685692  3.8382433535094176  0.6051784164185388 
C  0.9510334659475994  5.128478607796437  0.21585468295228194 
C  0.6208617880946186  5.420730736327316  -1.0824514159754552 
C  0.669330556566359  4.424726816933893  -2.0538897933498443 
C  1.005078645507899  3.133540941812697  -1.7315631939826335 
C  1.7228713490946366  3.5784023070330315  2.049681000778159 
C  0.5866957324224179  3.856133384190863  3.0370286988360866 
C  2.945961782999803  4.4174243944756135  2.4340393586594753 
C  1.0839254130713862  2.0598840640483753  -2.8026320648187766 
C  -0.04828351949222931  2.120567356445151  -3.822168911159828 
C  2.4424699855812984  2.124646991497328  -3.511435611460673 
C  2.0120037123886245  -2.087647377319398  -1.6241735024576553 
C  3.0413370390157146  -3.045934604483947  -1.7594697808548028 
C  3.396695318442969  -3.5200598006441126  -3.0259605710339517 
C  2.7287440519988637  -3.0510736678281756  -4.132991759791556 
C  1.708008954709478  -2.108180195925982  -4.027119328035877 
C  1.3670079623190774  -1.6472333694279295  -2.765774301355954 
C  1.4091963994062473  -2.9273452736059675  1.070904043066237 
C  0.16584332929612056  -3.5794281674296413  1.1103291953208803 
C  -0.025585499944217327  -4.67719744908087  1.8925825478855194 
C  1.0160286853916105  -5.166583070268819  2.6679696142476637 
C  2.2464323792668717  -4.538617576975707  2.6271899706580997 
C  2.4479839118595095  -3.4148829620267813  1.8325562972744645 
H  4.0659268051997275  2.3332113962202623  -0.6727865500377764 
H  4.4146310467380125  1.9129667445958367  1.8381384577460944 
H  5.737949710026472  2.0131949141614927  0.9438070374001728 
H  4.880724004232601  -0.18431740182914197  2.243935604271818 
H  6.147327445502299  -0.11137075312338503  1.2606084974521343 
H  4.660788824755782  -1.5627904600608318  0.15117231030902104 
H  2.8275983340629516  -0.4698246344236454  1.595780142498341 
H  5.5355163604665485  0.38591926500384266  -1.0724522913183803 
H  4.004874112279287  0.022153951501409114  -1.5944149382164068 
H  0.9315984783817401  5.81658669525096  0.8681900577695614 
H  0.35796741056116566  6.302881409729116  -1.3228390106596453 
H  0.462762634326102  4.645154403085465  -2.9517800999935764 
H  1.9635544286292124  2.6681079888059878  2.1273314961359917 
H  0.9170472785006287  3.8094059664938924  3.9199264757271335 
H  -0.0906442510265697  3.1882712670691533  2.9180160972552223 
H  0.21762321577959898  4.70607336311803  2.865197387137016 
H  3.2482547039023513  4.152242485054341  3.291721890187237 
H  2.6931951946651105  5.33948857797408  2.459457331256539 
H  3.624998432039018  4.293651143012919  1.79333244322624 
H  1.0304799678232583  1.2126071734033566  -2.36498341011828 
H  -0.8845740196746843  2.06227466642291  -3.367748520338122 
H  0.031643179189961845  1.4027958259323676  -4.4306832518479355 
H  -0.00508503503534552  2.949532022189531  -4.287174689498034 
H  3.1418576183768545  2.0860581342663105  -2.8637907634226223 
H  2.5110512377710696  2.947495990516313  -3.9878678381850494 
H  2.528129142920986  1.3999108286905364  -4.112443767049204 
H  3.4901488321885745  -3.3686877261015593  -0.9879303024124655 
H  4.089937701730053  -4.15869609542453  -3.119352920846888 
H  2.9678765195038608  -3.3757738092377565  -4.989942925934404 
H  1.260694167589828  -1.7913269883483418  -4.80638473601432 
H  0.6650421056133315  -1.0100464662491937  -2.682076293990991 
H  -0.5521954378799792  -3.242989914163244  0.5817086408857468 
H  -0.8815991451528631  -5.111880119116085  1.9172699047243809 
H  0.8843559989098151  -5.925842442445999  3.221632471405146 
H  2.9668716115552876  -4.8762573814972825  3.1509480111522983 
H  3.297813578724611  -2.994703972780234  1.810401065121912 
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

