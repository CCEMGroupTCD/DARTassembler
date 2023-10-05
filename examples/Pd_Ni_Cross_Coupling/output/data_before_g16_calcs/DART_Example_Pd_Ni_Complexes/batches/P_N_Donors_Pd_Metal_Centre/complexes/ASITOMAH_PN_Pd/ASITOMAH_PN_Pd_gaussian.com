%chk=ASITOMAH_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
C  3.1051791979997505  0.45496340973927113  -0.49575201846159594 
H  3.239130185186376  0.5106997144112844  -1.4540401470394677 
H  3.8887397500023324  0.8146805443769846  -0.0514255569566882 
C  2.873555798886171  -0.981260855599887  -0.07644641316389193 
C  3.960981991005583  -1.8446576274904094  0.203988393844389 
H  4.830239061129007  -1.5392730983074017  0.07813715079879072 
C  3.753682419107845  -3.1222772722034957  0.6569146140328782 
C  2.4126569517094216  -3.5497915564970555  0.9070864999711385 
C  0.512151394366916  1.355837004193256  2.5574207421929893 
H  -0.2513608447764546  1.0012871445495826  2.1611506478066556 
C  1.3337918063787604  -2.647716636464134  0.5775496631254265 
C  0.016809427047655534  -3.027537728871857  0.8660661622959244 
H  -0.680665813151627  -2.453974230378006  0.6413515377062512 
C  -0.2723835455449266  -4.2309127696241084  1.4751788086082516 
H  -1.1488234199395435  -4.4485191652289195  1.6992586451912277 
C  0.7803912994175324  -5.126042625974584  1.7513164664607983 
H  0.5813372333598126  -5.959348015698993  2.1151557939697403 
C  2.07807673454964  -4.810352391551488  1.5039912406226035 
H  2.752388567781706  -5.412346793350034  1.7237661977183594 
C  1.8016152094641558  3.061091819176044  -0.7117567668348227 
C  1.5491840859246482  4.207163423881094  0.028646079430226512 
H  1.262632923955304  4.133390914437118  0.9097216314825941 
C  1.7243042741677133  5.471014487472493  -0.5413836993783969 
H  1.552041113597237  6.233104342952265  -0.03679129388154731 
C  2.1486449683466384  5.595889290393137  -1.8454714030468773 
H  2.276008174671343  6.438305430463574  -2.2208200684548394 
C  2.3823632248581754  4.458009647652648  -2.5900468224706046 
H  2.6636365798680126  4.536043713427972  -3.4724277064091047 
C  2.202781790159266  3.206772190351605  -2.0399512909869144 
H  2.3520176908684105  2.4507854081544096  -2.560967088331566 
C  1.6552664804764148  1.6263593041278992  1.8095843127867866 
C  2.786379562161232  2.1403810648104167  2.431967174043715 
H  3.550393954422251  2.3186308932222746  1.9329072807551484 
C  2.7791956219095413  2.3845790274344996  3.7832394767469024 
H  3.5494988980083795  2.7032291784140443  4.1949919073555035 
C  1.6485548435521749  2.1667927643755327  4.539181586923525 
H  1.6385466056752682  2.376121370436015  5.445470023272229 
C  0.5518103958565757  1.640875703721655  3.941816496076987 
H  -0.19907337409893788  1.4615861938141632  4.459958470495702 
C  4.928080197003549  -4.048887894633733  0.9131956973634954 
H  4.602576907065133  -4.9024175713385585  1.2069897395485563 
H  5.492903097643524  -3.670655433324275  1.591439072653239 
H  5.430356721854814  -4.159248746302259  0.10266809450080662 
N  1.5904225617740713  -1.400948277060934  1.6984585296649457e-16 
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


