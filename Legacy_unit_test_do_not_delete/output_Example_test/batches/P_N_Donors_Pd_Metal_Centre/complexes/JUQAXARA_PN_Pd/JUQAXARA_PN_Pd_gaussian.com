%chk=JUQAXARA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5336702554656259  1.4628586902021679  2.326008497041909e-15 
N  1.5336702554656263  -1.4628586902021679  4.440892098500626e-16 
N  2.869529722474037  0.4451113746400297  -0.2432374890310056 
H  3.654766626905945  0.8041556647808725  -0.41514968067392033 
C  2.759316045574117  -0.9394775016239396  -0.17895756097970197 
C  3.909052531489319  -1.7462404158442448  -0.32147251877063493 
H  4.7657217859275125  -1.353250666492207  -0.44399437561491045 
C  3.7636813548904753  -3.114331328113713  -0.27906100170313297 
H  4.521449423830795  -3.67853730407571  -0.38212292669606185 
C  2.493220839635008  -3.6710997094038023  -0.08387304415335725 
H  2.3725260285844145  -4.612051822016792  -0.041280388042077454 
C  1.4199540200411573  -2.810360901593194  0.04483021411037158 
H  0.5562431583406056  -3.1850773172466056  0.17155576488988974 
C  1.6826789767863815  2.703062988232671  -1.3503376281941737 
H  2.5943451956476187  3.112273439665455  -1.3123552398352212 
C  1.5179448468158327  1.9946049417095482  -2.6959381494260506 
H  0.6509226223491846  1.5395286452369212  -2.7209444820176696 
H  2.2357770405108544  1.3375236655494016  -2.8086252747543616 
H  1.5621530039237514  2.6534862379383695  -3.4200373251782263 
C  0.6355677404755091  3.813285185959473  -1.1788519714580432 
H  0.7042773737818052  4.4449635026616106  -1.9249748973216587 
H  0.795405444813051  4.28559592904602  -0.3349007760461582 
H  -0.260916343813427  3.417482666633439  -1.1671925412921904 
C  1.9085819667216004  2.3564882928122035  1.5734088159477835 
H  1.0920509016805355  2.87628413850316  1.8238542563641604 
C  2.1512842514155786  1.3210373465621148  2.6808121882174007 
H  2.9732426803410625  0.8232396543603058  2.4883747210834675 
H  1.3948045260294735  0.6995979648304689  2.7188040282987394 
H  2.2406949650743417  1.7796236300446913  3.5426577456606636 
C  3.0727399792446617  3.3486021516015505  1.4813721217628857 
H  3.236822704508623  3.743227188774635  2.3638985798966523 
H  2.8487508953726306  4.057293053835399  0.8432170736772803 
H  3.8779575914046913  2.879250504994376  1.178537383176901 
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
