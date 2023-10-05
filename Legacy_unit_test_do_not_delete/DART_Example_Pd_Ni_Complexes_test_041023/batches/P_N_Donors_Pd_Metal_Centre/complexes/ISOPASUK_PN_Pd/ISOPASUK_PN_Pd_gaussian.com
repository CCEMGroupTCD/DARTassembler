%chk=ISOPASUK_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
C  2.2355770410181552  0.454963409739271  -1.4574064047854929 
H  1.6130608822276487  0.5106997144112844  -2.198171304322629 
H  3.090080298988578  0.8146805443769843  -1.7423929398463576 
C  2.392195526393042  -0.9812608555998871  -1.0047064603522302 
C  3.328229349408364  -1.8446576274904096  -1.6251740956567782 
H  3.816350159772541  -1.5392730983074019  -2.3553689078987583 
C  3.5261086380226088  -3.1222772722034957  -1.1680536951267448 
C  2.81470139639457  -3.5497915564970555  -0.004079892464395311 
C  2.7694523128369255  1.355837004193256  2.5125601304584126 
H  1.9640768351618698  1.0012871445495828  2.814803851627858 
C  1.847906112413948  -2.647716636464134  0.5771699742198514 
C  1.1810824377916025  -3.027537728871857  1.7489338346991592 
H  0.5473848969238555  -2.4539742303780057  2.1168955173424675 
C  1.4402334799258556  -4.2309127696241084  2.3714220102135783 
H  1.0203045063855631  -4.4485191652289195  3.172682444637556 
C  2.3498876218475058  -5.126042625974584  1.7738299902243506 
H  2.4870797674715526  -5.959348015698993  2.165212019238638 
C  3.0344102021288286  -4.810352391551488  0.643968895675591 
H  3.6489375101095822  -5.412346793350034  0.2899156963549806 
C  1.2027996679096518  3.061091819176044  -0.6332049602819622 
C  1.5841168211224192  4.207163423881094  0.04981412843702453 
H  2.0471434057233924  4.133390914437118  0.8523177600818438 
C  1.2776804183775268  5.471014487472493  -0.4617499048292565 
H  1.5373991004118854  6.233104342952265  0.003904793823256281 
C  0.5924937354262811  5.595889290393137  -1.6497014918486317 
H  0.3987779364543662  6.438305430463574  -1.9955080796654634 
C  0.19555440400057722  4.458009647652648  -2.321606206673678 
H  -0.271973773408166  4.536043713427972  -3.1210611009031384 
C  0.48419161745004313  3.206772190351605  -1.820065436188397 
H  0.19685973239340715  2.4507854081544096  -2.2795969399756264 
C  2.978594830726622  1.6263593041278992  1.1626598258251766 
C  4.197977815743722  2.1403810648104162  0.7385344134815255 
H  4.338329452003533  2.318630893222274  -0.16317517772307327 
C  5.19736184101186  2.384579027434499  1.6480507767698223 
H  6.018787027412611  2.703229178414044  1.351120036651553 
C  5.002589958953762  2.1667927643755327  2.9941046232345103 
H  5.669396702537143  2.3761213704360147  3.6079675243749096 
C  3.8247959054630574  1.640875703721655  3.4094293189547997 
H  3.707411120319608  1.4615861938141632  4.314149360654969 
C  4.502388096197347  -4.048887894633733  -1.869315810134767 
H  4.502915404634966  -4.9024175713385585  -1.4308331389492552 
H  5.384361440948837  -3.6706554333242756  -1.8352276250899915 
H  4.23613729785689  -4.159248746302259  -2.7849288369108556 
N  1.5904225617740713  -1.400948277060934  1.2328916564124948e-16 
P  1.5904225617740713  1.400948277060934  -1.160555310961094e-16 
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


