%chk=KEJOTUPA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.450537033308698  1.5453291930847604  2.537494656147078e-15 
N  1.4505370333086982  -1.5453291930847608  0.0 
C  0.9992871115304239  -2.2376681780645886  1.2386122629670027 
C  1.7372493256287973  -2.636799720884132  -0.9757790827733177 
C  2.7096929842627473  -0.851865059490788  0.3114719515394537 
C  2.483294665992427  0.4171414756704468  1.0259720251590838 
C  2.723094660773654  2.2278831141004387  -1.15703674517577 
C  0.9223197937010346  2.8770587230834335  1.1849975386939873 
C  2.170699492384158  3.0200360647325275  -2.317531178917209 
C  3.8583706696868605  3.0043425347658874  -0.4861253936202843 
C  0.5606556886703323  4.1875247513330685  0.49752228536774146 
C  1.905917629145047  3.114752258961555  2.3349423682887167 
H  0.7986247625027623  -1.5926325769379317  1.9045573845512695 
H  0.22086726449949312  -2.7491491968255195  1.045550590966484 
H  1.6884249668971743  -2.8120880259340586  1.5474995962685598 
H  2.0328268471646718  -2.2553217628943325  -1.7917982959850272 
H  2.413328327720386  -3.2016622006108832  -0.6220152181294005 
H  0.9457452547960772  -3.138625486209892  -1.1240630910495253 
H  3.166353583511439  -0.6734988471436416  -0.49755220223272434 
H  3.2402079906129133  -1.4191450500085843  0.8592972981558784 
H  2.040176902173853  0.24251459926923413  1.8436181842198376 
H  3.322019798197184  0.8298277743916662  1.195056919386536 
H  3.0919991679348695  1.4318039482469258  -1.513713119905188 
H  0.1315450346377136  2.5299876734973137  1.5716683379115421 
H  1.47015718775509  2.5386517691863832  -2.7252209377124887 
H  1.8454407684918877  3.855918571696305  -1.9996051543707183 
H  2.865812425308569  3.1761107872442707  -2.9472415182591343 
H  4.197018731675593  2.4908859369415057  0.2415999476053059 
H  4.553237149286613  3.150752824568515  -1.1175527073853662 
H  3.5303888390619926  3.8311462576991926  -0.16859983162000924 
H  -0.048016666005994635  4.0067741279384315  -0.2080726386469058 
H  0.15281659361692435  4.764909300847018  1.130839334347633 
H  1.3448792458091274  4.5871960496970825  0.15358537880984946 
H  2.114434940554746  2.291738984061079  2.7414345209559943 
H  2.692888278482474  3.520119806409188  1.9895680227479073 
H  1.501629811806857  3.6983444081526  2.9686390655773573 
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


