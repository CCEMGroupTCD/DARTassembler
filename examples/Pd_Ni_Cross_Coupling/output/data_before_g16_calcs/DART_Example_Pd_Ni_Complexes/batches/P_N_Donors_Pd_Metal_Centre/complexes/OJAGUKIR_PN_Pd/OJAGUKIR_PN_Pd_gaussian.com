%chk=OJAGUKIR_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5898296300232928  -1.4016211141032373  1.225976420596599e-16 
N  1.5898296300232928  1.4016211141032373  -1.7164908110038768e-16 
N  1.6316691832469596  2.747635490611154  0.2342106957849076 
H  0.881608625888068  3.1636696771203128  0.47879834050587367 
C  2.9070543248889975  3.1973412446538796  0.19922905403371471 
C  3.7238221369223368  2.113358005736938  -0.08184661010925424 
H  4.658101406260683  2.112108099480266  -0.17871432722375905 
C  2.8505230611758874  1.008266451810821  -0.19284242868611817 
C  3.249685161131457  4.672459252078327  0.364734534994106 
C  2.1626896105860918  5.394602272624638  1.1934691493133474 
H  1.3096700863482837  5.329364756522792  0.7376999996750784 
H  2.402867014153369  6.32968744422506  1.2959181447632522 
H  2.0937470738553303  4.981850271260747  2.0694166888379684 
C  3.3502460380574837  5.331808105845616  -1.0059003431071594 
H  2.494841087529983  5.2658044554517005  -1.4603645867476827 
H  4.030482533769643  4.883094920326583  -1.5321838760460234 
H  3.588742764658443  6.264086869950874  -0.8988297243595409 
C  4.604979857403775  4.768684432794366  1.101252729509528 
H  4.842001025969551  5.701351554348383  1.2208665498810045 
H  5.290953875099877  4.326073477695478  0.5777751946789083 
H  4.5332014732606725  4.339412007496352  1.9684091620224708 
C  3.0892785887254854  -0.44548409781875603  -0.4931969650530948 
H  3.8684312209768006  -0.7671834504205381  0.00600088369172001 
H  3.262033889759567  -0.5668464176984997  -1.448702419395473 
C  1.9606001535329416  -2.1554255369739534  1.6137410204109572 
C  3.2207606983185766  -2.7295704338385183  1.8334794530301868 
H  3.869566827157671  -2.7269442676349964  1.1538835212943355 
C  3.502093978185738  -3.3089031413937424  3.08183232365448 
H  4.34561110302124  -3.692357217956985  3.2368171334789007 
C  2.556373457254339  -3.322000157994941  4.075198386314898 
H  2.7550686428421507  -3.7079946928633363  4.909464664115917 
C  1.3238809719533458  -2.7712171844565314  3.852097064631654 
H  0.681276036267969  -2.7901483454764455  4.538146084582296 
C  0.9980145820166344  -2.1807943186241565  2.6279191913583597 
H  0.14564773091059835  -1.8118625047431678  2.487137417530668 
C  1.4563624779502742  -2.8089469992173735  -1.1463558881975848 
C  0.524306798567717  -3.801341831175897  -0.8260093133247473 
H  -0.02660440917428475  -3.7033849806272787  -0.07116564963347731 
C  0.41023931776216016  -4.944061603219821  -1.630669858152545 
H  -0.21429667125768503  -5.609508203313847  -1.4098616362895249 
C  1.2127078560982196  -5.091341497852455  -2.7463031929081 
H  1.1173290523622952  -5.848968540792239  -3.294766580568782 
C  2.149069024946263  -4.1382944428680855  -3.055197029782764 
H  2.719756173125095  -4.257654788054551  -3.793260866385335 
C  2.252631241661729  -2.9874416526459284  -2.269387305672109 
H  2.8721692053574674  -2.3222301588846217  -2.50591508439504 
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


