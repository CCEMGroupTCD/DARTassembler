%chk=YESOYAMI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4906405049508082  -1.3631180744895142  -9.464991743555295e-16 
N  1.490640504950808  1.3631180744895142  -4.444895748346378e-16 
N  3.821163363639565  1.7002164804160023  -0.22101058764468356 
N  2.9107686441832605  -0.44970840125465283  -0.1794232342847501 
C  2.731677813544815  0.9280171789034046  -0.13496336294692962 
C  1.3023511201780318  2.832773233122957  0.0043934364680217695 
H  0.49740037447449936  3.0522235548509746  0.5004884508400921 
H  1.1888898451502596  3.1455112438485857  -0.9071400438367104 
C  2.484477819653527  3.536479953730585  0.6328802444113972 
H  2.35447891324086  4.496852715635497  0.5953656802710331 
H  2.5668554130200674  3.2767305221689806  1.5641599955681538 
C  3.7242058479152256  3.162711481626065  -0.10793144907288123 
H  3.7073467111568954  3.5584456004051868  -0.9938178451428725 
H  4.50046806074009  3.5061525185848934  0.3615285823173103 
C  5.1593285464785215  1.1985520921302761  -0.5849773511870763 
H  5.7368084928914955  1.2382267182775761  0.19420879828620163 
H  5.540757738273691  1.772596585886285  -1.268556981394771 
C  5.113751486735381  -0.20741141966434654  -1.0919218798362273 
H  6.0108715325398645  -0.5769904931343127  -1.1185358210868062 
H  4.7543709302123025  -0.21978699835830343  -1.992677096665506 
C  4.245894520187006  -1.0392241813171679  -0.18415824035455278 
H  4.205937260425451  -1.9531988899766302  -0.5050465532037376 
H  4.611897066976063  -1.0444695156237191  0.7145202795912414 
C  1.7061180846004222  -2.2779776624727655  1.5498514338295495 
C  2.382217248984893  -3.491402061250132  1.6592789823281668 
H  2.723880433399205  -3.8992155077912054  0.8968112051492352 
C  2.5538247459115317  -4.102803458236311  2.907469808589079 
H  2.9904321770909945  -4.9212364178445975  2.971930895034263 
C  2.0666486240266373  -3.479807188952719  4.043534210251213 
H  2.182568489815684  -3.875745248245287  4.876848081442075 
C  1.4104454494294252  -2.271870946440026  3.9392144725556566 
H  1.0963511308288871  -1.8529306978384048  4.706389893906276 
C  1.2162122498270573  -1.6753099623121217  2.708023932310713 
H  0.7567720685715934  -0.8682819318078848  2.6511214909887455 
C  1.5966056760122724  -2.5379197112440712  -1.3745589244556873 
C  1.1453253018274823  -3.8572039894409365  -1.263734255819972 
H  0.8014531128497274  -4.162614326073666  -0.4549928855874667 
C  1.207028600716219  -4.713868699456842  -2.3558710272545755 
H  0.9383836413302717  -5.599839762860338  -2.272385216434991 
C  1.6717304029889577  -4.238843720243289  -3.5706910814550676 
H  1.7102639437516367  -4.8054877288991555  -4.307392957175988 
C  2.0795977226893174  -2.917958891044401  -3.6914613489577386 
H  2.365834147509371  -2.597774911264477  -4.51464322327139 
C  2.064373817532108  -2.0786500252820166  -2.603912673925456 
H  2.3666375391360917  -1.2038759266064043  -2.687709624277927 
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


