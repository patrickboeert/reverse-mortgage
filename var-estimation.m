
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% real estate, 1-year, cpi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model:  3-D VAR(1) with Additive Constant
          Conditional mean is AR-stable
  Series: Index Wohnen
  Series: WZ9808: Zinsstrukturkurve (Svensson-Methode) / Börsennotierte Bundeswertpapiere / 1,0 Jahr(e) RLZ / Monatsendstand
  Series: UJFB99: Verbraucherpreisindex / bis 1994: Westdeutschland

pvalue =
  1.0e-003 *
       NaN       NaN       NaN       NaN       NaN
         0       NaN       NaN       NaN       NaN
         0    0.0012       NaN       NaN       NaN
         0    0.0000    0.1595       NaN       NaN
         0    0.0000    0.0009    0.4942       NaN
         
reject =
   NaN   NaN   NaN   NaN   NaN
     1   NaN   NaN   NaN   NaN
     1     1   NaN   NaN   NaN
     1     1     1   NaN   NaN
     1     1     1     1   NaN         
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% real estate, 1-year, 10-year, alpha=0.01  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

Model - Information
  Model:  3-D Time Series with Additive Constant
  Series: Index Wohnen
  Series: WZ9808: Zinsstrukturkurve (Svensson-Methode) / Börsennotierte Bundeswertpapiere / 1,0 Jahr(e) RLZ / Monatsendstand
  Series: WZ9826: Zinsstrukturkurve (Svensson-Methode) / Börsennotierte Bundeswertpapiere / 10,0 Jahr(e) RLZ / Monatsendstand
  
pvalue =
       NaN       NaN       NaN       NaN       NaN
         0       NaN       NaN       NaN       NaN
         0    0.0000       NaN       NaN       NaN
         0    0.0000    0.0003       NaN       NaN
         0    0.0000    0.0000    0.0191       NaN
         
reject =
   NaN   NaN   NaN   NaN   NaN
     1   NaN   NaN   NaN   NaN
     1     1   NaN   NaN   NaN
     1     1     1   NaN   NaN
     1     1     1     0   NaN  
     
pvalue =
       NaN       NaN       NaN       NaN       NaN
         0       NaN       NaN       NaN       NaN
         0    0.0000       NaN       NaN       NaN
         0    0.0000    0.0003       NaN       NaN
         0    0.0000    0.0000    0.0191       NaN   
         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With suggested p=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model - Parameters and Standard Errors
  Model:  3-D VAR(3) with Additive Constant
          Conditional mean is AR-stable
          Standard errors without DoF adjustment (maximum likelihood)
  Series: Index Wohnen
  Series: WZ9808: Zinsstrukturkurve (Svensson-Methode) / Börsennotierte Bundeswertpapiere / 1,0 Jahr(e) RLZ / Monatsendstand
  Series: WZ9826: Zinsstrukturkurve (Svensson-Methode) / Börsennotierte Bundeswertpapiere / 10,0 Jahr(e) RLZ / Monatsendstand

       Parameter          Value     Std. Error    t-Statistic
  -------------- -------------- -------------- --------------
            a(1)       -2.06795       0.973969       -2.12322 
            a(2)      -0.939542       0.804425       -1.16797 
            a(3)     -0.0999303       0.554255      -0.180297 
      AR(1)(1,1)        1.16709       0.174603        6.68427 
           (1,2)      -0.275184       0.276801      -0.994158 
           (1,3)      0.0199303       0.389339      0.0511901 
           (2,1)       0.350109       0.144209         2.4278 
           (2,2)       0.165653       0.228617       0.724586 
           (2,3)       0.703528       0.321564        2.18783 
           (3,1)      0.0152082      0.0993608        0.15306 
           (3,2)      -0.193968       0.157519        -1.2314 
           (3,3)        1.00756        0.22156        4.54758 
      AR(2)(1,1)      0.0345545       0.239547        0.14425 
           (1,2)      -0.457591       0.290161       -1.57703 
           (1,3)        1.12841       0.463649        2.43375 
           (2,1)       0.140179       0.197847        0.70852 
           (2,2)      -0.401981       0.239651       -1.67736 
           (2,3)        0.48349       0.382939        1.26258 
           (3,1)       0.348837       0.136318        2.55899 
           (3,2)      -0.113854       0.165121      -0.689518 
           (3,3)      0.0348604       0.263848       0.132123 
      AR(3)(1,1)      -0.421124       0.166446       -2.53008 
           (1,2)     -0.0621675       0.260968      -0.238219 
           (1,3)      -0.111817       0.376279      -0.297165 
           (2,1)      -0.154032       0.137472       -1.12046 
           (2,2)      -0.115822        0.21554      -0.537359 
           (2,3)     -0.0460021       0.310778      -0.148023 
           (3,1)      -0.222057      0.0947194       -2.34437 
           (3,2)     -0.0532228       0.148509      -0.358382 
           (3,3)       0.202068       0.214128       0.943677 
          Q(1,1)       0.795293                
          Q(2,1)        0.20091                
          Q(2,2)        0.54251                
          Q(3,1)      0.0456621                
          Q(3,2)       0.256017                
          Q(3,3)       0.257546           