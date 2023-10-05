%chk=VERUSOLI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.4628804069369437  1.5336495411273074  -9.90144842482718e-17 
N  1.4628804069369437  -1.5336495411273074  1.1102230246251565e-16 
C  1.2691915704031649  -2.4069167798835363  -1.152745602601585 
C  1.4671427349562094  -2.1080727234616816  -2.469172515268381 
H  1.7959898318717349  -1.2978693903803489  -2.761503536397048 
C  1.1816754739162287  -3.0986081354068142  -3.387603721487456 
H  1.3153218353898608  -2.956812556189834  -4.224973050863323 
C  0.7103938744536706  -4.336023610134536  -2.9968888677739294 
H  0.4783026470191786  -4.9214021343181615  -3.6107272386603952 
C  0.5157269413291656  -4.628380486967313  -1.6654486168739893 
H  0.26524745407145334  -5.535921921744466  -1.458790588658918 
C  0.8099769346079304  -3.6563588842325476  -0.7100364360470931 
C  0.7457954867964456  -3.639538379085695  0.7293542414145457 
C  0.31488770238102215  -4.562551829707199  1.6756557691011436 
H  0.0025378949514864857  -5.333322115421303  1.364162375176115 
C  0.3348136105951496  -4.194684214016263  2.995776924718191 
H  0.04082361831692971  -4.80150381330746  3.63649511421835 
C  0.7758167503228137  -2.954112050051415  3.4073375724311754 
H  0.7784902833285084  -2.7329984573392654  4.310192467383431 
C  1.2134180678180775  -2.041689304518939  2.4635045626453382 
H  1.4623412036707322  -1.2460662744340554  2.711970911344354 
C  1.1742962531714074  -2.38153993924536  1.1518613757630174 
C  2.790287328149698  -0.9019489380756617  0.019765971281469855 
C  3.90402511063286  -1.7251906345273507  0.06356929413303788 
H  3.7401082587780516  -2.6584947086431887  0.15066797906449586 
C  5.1693191072486515  -1.196031052603912  0.00622551516204542 
H  5.858962397590094  -1.7292917018385796  0.18170876810911193 
C  5.344982480769123  0.16930982146621143  -0.11509264834549376 
H  6.199899716333026  0.5303359947554571  -0.17827431131752586 
C  4.23472773656912  0.9940346175042163  -0.1426257880409005 
H  4.389884975989475  1.9311582213384906  -0.24112229656847262 
C  2.9445781358501977  0.48025923801248593  -0.05386584666568746 
C  1.6302419667820383  2.6147779968377516  1.4701341641290515 
H  0.8181710465838018  3.1632183298109147  1.4972692288367293 
C  1.616774473420233  1.7749304107632748  2.754312418417519 
H  0.8179944105621813  1.200127676266304  2.7364614353831933 
H  2.450307920623493  1.2560909245172138  2.719927539395186 
C  1.6477132801986516  2.6585041501743922  3.9870841033852824 
H  0.838658465221378  3.0497706546190257  4.048089667572317 
H  1.6901026534833163  2.0777323922763706  4.782396511348718 
C  2.793357704113279  3.6353981339842405  3.979560251349724 
H  3.6639722898182  3.0952693185111104  3.961155473343406 
H  2.6662052417279494  4.133504392626929  4.704766943548794 
C  2.7741212318810438  4.473416772228864  2.7355377635743614 
H  3.455154534540363  5.131014879198391  2.6598318828314937 
H  2.004659379868391  4.959143018690908  2.6452778511261603 
C  2.8064140192134954  3.600489581825284  1.4870245253023882 
H  3.6420098930693294  3.1062602193316535  1.4605730065512441 
H  2.766723694471411  4.164352799470402  0.6981979791537516 
C  1.6833434907744562  2.622461264549187  -1.4651597556594131 
H  2.643130705084773  2.9970096568616653  -1.424169939087316 
C  0.7218085088934023  3.810247027299782  -1.4736932441601804 
H  0.8733963346974886  4.359340697279304  -0.6895137987164596 
H  -0.19328906786940347  3.4901201457443523  -1.4485054552845607 
C  0.9428275986458962  4.641530270838517  -2.7375661574614782 
H  0.3099046669739509  5.376615855757306  -2.753776197819275 
H  1.8383504539251638  5.0167664897839455  -2.7238922529710266 
C  0.7707081574936647  3.807305628447551  -3.9745432081803242 
H  0.9116339979311288  4.360148861052554  -4.759950572329585 
H  -0.1358283565809173  3.462975052368328  -4.008916643075285 
C  1.7580771160513946  2.6421333892366987  -3.988252790559912 
H  1.612751132390034  2.1003578937807714  -4.778832521514568 
H  2.6663655838638407  2.9834396009643434  -4.019362741942464 
C  1.573933238534027  1.791350159424686  -2.742688051879668 
H  2.2470127453213222  1.0940387846148236  -2.729635308728712 
H  0.7029122788618873  1.3656746465430851  -2.773073172330575 
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


