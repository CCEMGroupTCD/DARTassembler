%chk=MIMIBIGO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5978050538160153  -1.392522534826635  1.7053482650159972e-16 
N  1.5978050538160153  1.392522534826635  -1.7053482650159972e-16 
C  2.5929046909895224  1.0064872308687547  -0.7003623820532702 
C  2.746588439502177  -0.39486579660159926  -1.0058577958368022 
C  3.550095935024414  -0.8973306418740019  -1.962779670341189 
C  4.54279658890763  -0.21068785076180097  -2.774842586126349 
C  5.414507233156492  0.7517503909768514  -2.277606461591661 
C  6.368030042866021  1.3177095336761941  -3.1014049443279714 
C  6.446526502182868  0.9627639743488298  -4.416970440906768 
C  5.608806145690355  0.023584113849969786  -4.923387417272523 
C  4.671589589433539  -0.5758302648221246  -4.1144960900107375 
C  1.052753018796467  -2.79649786815465  -0.9917237783413092 
C  0.48011608369669156  -2.5488150469621447  -2.2302779921913523 
C  0.018359031122132885  -3.5987230688686025  -3.005130507147664 
C  0.08887810451209521  -4.888670272648415  -2.516406820644912 
C  0.634659325908216  -5.134377283804682  -1.3016477074705846 
C  1.1303139211597149  -4.099562611897489  -0.5251286433055776 
C  2.5851633720625724  -2.015500560893553  1.390420317899284 
C  3.960496243265212  -2.0734494399830004  1.3444826340022018 
C  4.670392074429287  -2.4993417075532585  2.4504594037942926 
C  4.028568813725572  -2.875369510984956  3.5807123325964993 
C  2.6696895639020517  -2.8383057804461713  3.640244460709162 
C  1.9316443509338217  -2.3944298541747497  2.5540442003099155 
H  1.5375493376749685  2.2275207132545045  0.19447398462126114 
H  3.214864613617493  1.6238392007649007  -1.0108984667680163 
H  3.4588308379709014  -1.809244257805606  -2.1257984702949395 
H  5.35640629067091  1.0142418270996478  -1.3872851395710792 
H  6.963101216704217  1.943618190952257  -2.7559138434096844 
H  7.076843852099502  1.3667493002093105  -4.968965753281322 
H  5.669312129563577  -0.2157317045051615  -5.820757605656639 
H  4.11533735500124  -1.2326527557514155  -4.465987593604342 
H  0.4063510394481087  -1.6737514917775924  -2.541895406073791 
H  -0.33696098947892605  -3.437791655268248  -3.848556119206407 
H  -0.24220473845087165  -5.593097608987446  -3.0258868486282067 
H  0.6769430651605896  -6.007474587146609  -0.9846981936803083 
H  1.513847411304946  -4.278009971839809  0.3034749661885227 
H  4.411801981441245  -1.825151540773028  0.5697929839246318 
H  5.59889215150192  -2.5304899359770077  2.4176898334257984 
H  4.517335460558682  -3.1591606598651416  4.319491953837362 
H  2.233494763136668  -3.112171093546816  4.413727082405042 
H  1.004351843158262  -2.3510082732036905  2.605752529709591 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Pd 0
lanl2dz


