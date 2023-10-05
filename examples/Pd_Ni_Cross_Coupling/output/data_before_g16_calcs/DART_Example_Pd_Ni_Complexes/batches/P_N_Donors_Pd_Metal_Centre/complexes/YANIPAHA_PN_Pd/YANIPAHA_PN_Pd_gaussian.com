%chk=YANIPAHA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5832911782107546  -1.4090028548587124  1.7255308361922047e-16 
N  1.5832911782107546  1.409002854858712  -1.7255308361922042e-16 
N  2.9978428436500044  -0.4205920221924786  -0.19533487059696208 
C  2.7680701680924114  0.9251119225948674  -0.19718704523335084 
H  3.5912326512464343  1.516530312647427  -0.36931209929390074 
C  1.397197928371946  2.8329404168597385  -0.11906100069714323 
C  0.9595928484661151  3.3630076253288737  -1.3363741356572492 
C  0.7705076877426776  4.741865438356384  -1.4124082296361125 
H  0.4737826015436353  5.118380406491914  -2.2325874490693383 
C  1.0000556286430675  5.584231879210058  -0.3360098033310902 
C  1.4025636709444078  5.014377788970819  0.8664540222039594 
H  1.5486981482596067  5.576903966449344  1.617272118810066 
C  1.5966371429635589  3.6428396310912303  1.0026342665079564 
C  0.7052826379005929  2.47772764473928  -2.5226902498399113 
H  1.5496838238826554  2.0783193695470827  -2.8185107643755156 
H  0.3223599522630125  3.0102327084534606  -3.250718322061627 
H  0.07703203381646895  1.768849116526394  -2.271444810664043 
C  0.8380806492816946  7.0773844854984  -0.45912625387777956 
H  1.7000926836569281  7.483319500482706  -0.6873501487748471 
H  0.5202537029174252  7.441983789345457  0.3932281910563813 
H  0.18632100033539367  7.2791882426870345  -1.1627151562845486 
C  2.0104246437084066  3.054858100647083  2.319218465898163 
H  1.429525865109406  2.294398076520698  2.5351604505326826 
H  1.933228378757235  3.7357171505006033  3.0204868816350747 
H  2.9406364980048645  2.7483120468115056  2.262349906885296 
C  4.346754928044385  -0.9184758040321861  0.05375179776530967 
C  5.0661419235078595  -1.5627693353543994  -0.9602948117425794 
C  6.323778634755594  -2.078719457451844  -0.6310583860092991 
H  6.823095420117275  -2.528052996994946  -1.302782524400776 
C  6.86448253700624  -1.9590938375187306  0.632206044403663 
C  6.124014484570338  -1.2995188715807153  1.6045868790716828 
H  6.48940518154771  -1.2025017580038184  2.475894774270382 
C  4.867667592152769  -0.780462907656969  1.3421467864258814 
C  4.548615735095544  -1.717415395331227  -2.360117443521739 
H  3.9808259592532833  -2.515865456039452  -2.4130446671785983 
H  5.302585741302325  -1.8119688978600488  -2.978653590323302 
H  4.02247943324  -0.9267388603914067  -2.603329982189743 
C  8.235117618101217  -2.5098920224364045  0.9543840758757256 
H  8.315773198567017  -3.4160351083764455  0.5932821579317874 
H  8.357050961272183  -2.5333914389818175  1.9271659719064784 
H  8.920088265735263  -1.9349228006165986  0.5540554910360365 
C  4.1218672125594615  -0.07140001147129507  2.4448911552853705 
H  4.239659929366641  0.896069366109698  2.3478593599282247 
H  4.472325210265986  -0.35694677144598913  3.3136034492650888 
H  3.1679578644761857  -0.2919576802711665  2.3878167269226127 
C  2.048139876483148  -2.6583923035241828  1.1935274246023424 
C  1.595180660403151  -2.554238046849723  2.51348668889228 
H  0.9767404647275856  -1.8741264227351704  2.749972974154304 
C  2.051458692692661  -3.4410085283228002  3.470470252198373 
H  1.7562319339825812  -3.3596472063036935  4.370042142722355 
C  2.9324298891528073  -4.447954238267632  3.123719494115131 
H  3.232270954303311  -5.062401612887617  3.7842567371771922 
C  3.382610760701042  -4.568440523596715  1.812032813562299 
H  3.9861985386610623  -5.262733843604286  1.5760516912449352 
C  2.9468371090345604  -3.6728609187374097  0.8557973949161791 
H  3.2602184368923757  -3.746408221812527  -0.03796392532568536 
C  1.279314536779286  -2.2072833914249697  -1.585903854716164 
C  1.2799129992956768  -1.3972990992963048  -2.7242701614678513 
H  1.4743863559966148  -0.47114742366265083  -2.6470258914250686 
C  0.9974109952543508  -1.9390501201461134  -3.962239025826262 
H  1.0121140803854327  -1.3912431144055908  -4.737753696060949 
C  0.6954649359432619  -3.2867275030925165  -4.066719783677851 
H  0.5077966701689842  -3.663240210764211  -4.91911048223281 
C  0.6664786922122374  -4.085167180571248  -2.944853045111252 
H  0.45199411603678863  -5.007068850621628  -3.0291696943636786 
C  0.9474696246593949  -3.555236767424658  -1.6930061442863933 
H  0.91384290994633  -4.105204840468558  -0.9188756819675442 
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


