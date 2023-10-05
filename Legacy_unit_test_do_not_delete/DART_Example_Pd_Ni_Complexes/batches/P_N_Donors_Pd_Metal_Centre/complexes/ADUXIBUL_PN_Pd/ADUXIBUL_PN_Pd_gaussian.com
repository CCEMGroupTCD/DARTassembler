%chk=ADUXIBUL_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5832911782107546  -1.409002854858712  1.7255308361922042e-16 
N  1.5832911782107546  1.409002854858712  -1.7255308361922042e-16 
N  3.0061272869949436  -0.42059202219247815  -0.1210352568797772 
C  2.7767664422342753  0.9251119225948675  -0.13491026586900465 
H  3.607809139111424  1.516530312647427  -0.2637184328259155 
C  1.4036841340259854  2.832940416859738  -0.1286371999142732 
C  1.0303880236386203  3.363007625328873  -1.3671845298437568 
C  0.8455413145718297  4.741865438356383  -1.4530103744371665 
H  0.5921477435304985  5.118380406491913  -2.287594957882635 
C  1.0184403270579834  5.584231879210057  -0.3660735037750062 
C  1.3574646519194198  5.014377788970818  0.8558080301675812 
H  1.464104073961912  5.576903966449343  1.6132452442770058 
C  1.5441450296830619  3.64283963109123  1.0019586649191563 
C  0.838513324620398  2.47772764473928  -2.5651844070741334 
H  1.6972393376283759  2.0783193695470827  -2.8164069663546356 
H  0.4942174665465999  3.01023270845346  -3.3122553670451333 
H  0.19797454584364194  1.768849116526394  -2.347163387165677 
C  0.8631307459450518  7.077384485498399  -0.4974983429373896 
H  1.735905738562962  7.483319500482704  -0.6802952408063717 
H  0.5011305854312817  7.441983789345455  0.3370542026954555 
H  0.24908730895538445  7.279188242687034  -1.2342334656745124 
C  1.88846075598971  3.0548581006470825  2.338394495954281 
H  1.2970565487190433  2.294398076520698  2.5236386465454936 
H  1.7746687327161652  3.735717150500603  3.0346617073551068 
H  2.820374057685068  2.748312046811505  2.330287400228226 
C  4.340154545282418  -0.9184758040321856  0.19830665068484749 
C  5.111626744854108  -1.5627693353543988  -0.7767004368711494 
C  6.350309005526944  -2.078719457451843  -0.38209559832748374 
H  6.884096819702015  -2.528052996994945  -1.026766940702725 
C  6.824157739918131  -1.9590938375187295  0.907735827856393 
C  6.033813992398779  -1.2995188715807149  1.8400309447935985 
H  6.35310320224769  -1.2025017580038178  2.7292678143324487 
C  4.792923932938571  -0.7804629076569685  1.5121984011862502 
C  4.668070864101287  -1.7174153953312261  -2.2016898883350193 
H  4.103829221274892  -2.5158654560394513  -2.2842603979349914 
H  5.453379261332153  -1.811968897860048  -2.7799186115910204 
H  4.155384374543257  -0.9267388603914063  -2.4721049591432624 
C  8.176052917912745  -2.5098920224364036  1.3012058231627135 
H  8.275496576885423  -3.4160351083764446  0.9448199697773817 
H  8.246907684930042  -2.5333914389818166  2.279036053526328 
H  8.88103641637904  -1.9349228006165975  0.9372744685775742 
C  3.9904324658097394  -0.07140001147129484  2.5743993211812617 
H  4.113142003588356  0.8960693661096981  2.483665298999912 
H  4.294945284404995  -0.35694677144598863  3.4602626295867642 
H  3.0408174621206827  -0.29195768027116603  2.467479353437956 
C  1.9850384184290428  -2.658392303524182  1.2162200378798351 
C  1.4636186369139321  -2.554238046849722  2.510664290160084 
H  0.8336492561553851  -1.8741264227351697  2.714459820167382 
C  1.8691867061235337  -3.4410085283227994  3.490216087884742 
H  1.5272845903259404  -3.3596472063036926  4.373104171671401 
C  2.7670980946789623  -4.447954238267631  3.1900470096535103 
H  3.0319583898207942  -5.062401612887615  3.8653714782462 
C  3.285310385706449  -4.568440523596714  1.9037185965040233 
H  3.9004212652859738  -5.262733843604284  1.6996502216494636 
C  2.900179441712665  -3.672860918737409  0.9259870345001588 
H  3.2599071446999  -3.7464082218125263  0.04985169454751452 
C  1.3627309209477303  -2.2072833914249688  -1.5996393368055222 
C  1.4229060525106396  -1.3972990992963044  -2.7364142309891712 
H  1.6130702374980803  -0.4711474236626504  -2.6490978724471717 
C  1.2055814919628043  -1.9390501201461126  -3.9874715145661095 
H  1.260851728873877  -1.3912431144055903  -4.761153868877866 
C  0.9095173396031448  -3.2867275030925156  -4.107611720918082 
H  0.7667169489662451  -3.66324021076421  -4.9686560457908975 
C  0.8218568519717137  -4.085167180571247  -2.988799484510242 
H  0.6120790119194198  -5.007068850621626  -3.084225836118762 
C  1.0369460914872377  -3.555236767424657  -1.72396226721034 
H  0.9628506030157558  -4.1052048404685575  -0.9526526100553417 
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


