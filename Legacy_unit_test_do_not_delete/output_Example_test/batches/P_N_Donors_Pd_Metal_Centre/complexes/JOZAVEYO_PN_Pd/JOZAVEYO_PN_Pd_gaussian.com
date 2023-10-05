%chk=JOZAVEYO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.0326855125347283  -0.47796135277049734  -0.49692971958742105 
H  3.799030190040134  -0.7730962753971786  0.035633396583139465 
C  3.265003440798926  -0.8643134374329426  -1.9558536233206083 
C  4.254826555405462  -0.3812138788570034  -2.721519874564813 
H  4.79040648158276  0.2972159191131105  -2.3804990366352694 
C  4.530538519150822  -0.8920791178795549  -4.105914374194958 
H  5.489626507194416  -0.9331317084675589  -4.247704662989918 
H  4.157276054027458  -0.27853854167777437  -4.7558691599178236 
C  3.931948096418713  -2.2744274790997485  -4.3152250620378005 
H  4.425301315786155  -2.925814809443141  -3.790204874989819 
H  4.0075427589751165  -2.522757751303537  -5.249958726190453 
C  2.4614430377591203  -2.293985097628476  -3.898527675420011 
H  1.9521914333387969  -1.7051220854901055  -4.477746672286212 
H  2.1074760624242157  -3.1914775614493873  -3.997894660666301 
C  2.3078459713955306  -1.8465120398227632  -2.462913909902781 
C  1.373029800943864  -2.2845557143257196  -1.582486282729464 
C  0.32103353589214256  -3.28605312932188  -1.8118521173109645 
C  -0.3895015965107249  -3.359921269306674  -3.019210064862731 
H  -0.2503119547837216  -2.7288537456662336  -3.6890357916321967 
C  -1.3044963369630915  -4.3912934927666445  -3.1971467435432266 
H  -1.7733285690244172  -4.475697405044666  -3.9958123918420325 
C  -1.5058882930974895  -5.288088526236737  -2.1667004092065865 
H  -2.1060033470178308  -5.993479672058642  -2.2600013387532902 
C  -0.7982391805455986  -5.113462951332219  -0.9967496828536233 
H  -0.9688454151829862  -5.701291561519233  -0.2949437426190997 
C  2.8205881982469503  0.9926521644580755  -0.2726111911139185 
C  3.8773425213656454  1.8771105499886784  -0.33463920764323807 
H  4.743106754395669  1.5571739600286214  -0.45563330341114283 
C  3.6567163511588574  3.242748830472846  -0.21818651383099552 
H  4.366224335756403  3.8430379768998284  -0.2564253746180194 
C  2.3515159255173237  3.6942620767680276  -0.042597398777714957 
H  2.167016529037886  4.605025036770059  -0.013286250948578271 
C  1.3329831327967698  2.764254202556872  0.08798802942214207 
H  0.46713183389485824  3.065929689618198  0.24262026634642186 
C  2.011927517966126  -2.5792074167578045  1.3392287036163686 
H  1.2055077767799371  -3.075264552225777  1.5925613098409248 
C  3.082456803694625  -3.6111383944479707  0.9379939522786173 
H  2.7718818498283038  -4.124870903381496  0.17601181515478526 
H  3.897892620363595  -3.1532992758220804  0.6830822413877647 
C  3.3619750962298673  -4.550386776999513  2.1160086635672988 
H  4.056312008039926  -5.179437512417847  1.869847405105031 
H  2.5590575122056913  -5.053046981013205  2.325805142483675 
C  3.8067825645834428  -3.748975107779226  3.3529583604161455 
H  3.947094066766551  -4.3567573588258375  4.096868138324654 
H  4.647106028639557  -3.3063566572872993  3.165003996165447 
C  2.766075597811492  -2.717663686010706  3.733382887283374 
H  3.088146212400933  -2.2000328719143463  4.488158568870933 
H  1.9507747065271328  -3.1680214990678355  4.006150181674288 
C  2.457942092333099  -1.7745770441273356  2.5624939990011537 
H  1.7566164653856673  -1.1548007614243518  2.8183445498813016 
H  3.250109315025104  -1.2595407556756055  2.3422767588261033 
N  1.5532379461949806  1.4420651450263955  8.832590856395058e-16 
N  0.11305523835277453  -4.156975091018009  -0.8023166477072285 
P  1.55323794619498  -1.4420651450263955  4.440892098500626e-16 
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

-Pd 0
lanl2dz
