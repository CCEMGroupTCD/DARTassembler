%chk=LETUGEYO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4960737665636685  1.5012872093640177  5.327208233055968e-16 
N  1.4960737665636674  -1.5012872093640177  1.1102230246251565e-16 
C  1.2402199987950144  -2.7108104662558388  -0.8210258403706935 
H  0.3091370812406127  -3.017117085756585  -0.6845501814111717 
H  1.8506121014004655  -3.43939159861849  -0.5433647076169993 
C  1.4607669726206638  -2.3884584259656862  -2.2782176789620077 
C  2.7507941521411405  -2.2012351496686287  -2.765086851236262 
H  3.490249589222844  -2.2547467925434326  -2.1724845810354485 
C  2.9639093692948855  -1.9371713499324401  -4.099625357194969 
H  3.851371923872417  -1.82359485774862  -4.417737058718091 
C  1.8968111691028084  -1.8352137504502974  -4.99415877694781 
C  0.6171797871557085  -2.0016746038041786  -4.502063785470626 
H  -0.122833758675432  -1.9275061107233844  -5.092770520340462 
C  0.39264081126141503  -2.274346350054127  -3.1579304375094184 
H  -0.49616331697200566  -2.3837908513961494  -2.840284388012127 
C  2.1508482899893675  -1.550479342314571  -6.448806249030087 
H  1.3421293017931124  -1.7463820304648228  -6.966085269391777 
H  2.8872750932085403  -2.113342047551346  -6.76847487206636 
H  2.3908757174072335  -0.6066487540287626  -6.560724856523308 
C  2.660449456272941  -1.4277238736404128  0.5212202786263349 
H  3.2167656798559956  -2.1900649652459525  0.41015966712087576 
C  3.2402129626143448  -0.3069045364045866  1.271200458359747 
C  4.329661447985529  -0.6110315529755337  2.087506507286338 
H  4.647895431062286  -1.5052597004122341  2.122840093957218 
C  4.958240287185259  0.3647627892813272  2.848393295533986 
H  5.682556722007058  0.13431182837121458  3.4184002488798884 
C  4.524014401190639  1.670007580403016  2.770235217372652 
H  4.94337287997362  2.3425888284826217  3.2932207439453105 
C  3.4682653555780885  2.0044635405072175  1.9233454645482038 
H  3.1898183735152497  2.9105993594996984  1.8600622462440257 
C  2.816871562557692  1.0345388679211802  1.1704289247356123 
C  2.316987460952678  1.6093946580838097  -1.6199804924430294 
C  3.6941360623519257  1.52426523101426  -1.76770946831066 
H  4.247950108083164  1.456789211523684  -0.9984067453566082 
C  4.2683410172478355  1.5368950128985888  -3.0312325818240438 
H  5.211658820364797  1.4739973778484634  -3.125622304928391 
C  3.47205400820339  1.6389826241458247  -4.148369891583675 
H  3.865917393746189  1.64989153611007  -5.012708916022879 
C  2.1006696386036343  1.725546453631225  -4.0154721706012655 
H  1.5534459692286593  1.8028277588180917  -4.787923434788497 
C  1.5229999139399417  1.7013162698441993  -2.759174432244354 
H  0.5776662302944443  1.74866127788574  -2.673329416205001 
C  1.0975475347188184  3.204939444623846  0.49895139425173096 
C  1.5176222344970247  4.315546711907716  -0.21420365513271097 
H  2.008745347600981  4.205572868148705  -1.0202749804284088 
C  1.2202187673337308  5.598657477113019  0.25432343259383716 
H  1.5022638867704783  6.360752922850335  -0.23609191401387586 
C  0.5118036165715488  5.754689173472595  1.4352109055995372 
H  0.315196290539403  6.62717298724923  1.7553390893763927 
C  0.09267892066885608  4.654218490917536  2.146883678804507 
H  -0.3903393520019378  4.768318507177961  2.957370461272618 
C  0.3738442245134559  3.3785869508264974  1.684581801413847 
H  0.07526342926764329  2.6207808475679806  2.173101635785517 
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

