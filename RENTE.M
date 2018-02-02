classdef rente < handle

   properties
       
       % Tarife + Beispiele fuer Versicherung
       Versicherung = []; 
       
       % Zahlungen
       RenteGarantiert = {};           % Versicherungsleistung (Garantiert + Ueberschuesse), {tarif,beispiel}(:,szenario)
       RenteFlexibel = {};             % Versicherungsleistung (Flexibel), {tarif,beispiel}(:,szenario)
       Zinssatz = 0.0463;              % Zinssatz Hypothek
       Zinsen   = {};                  % Zinszahlungen fuer Hypothek {tarif,beispiel}      
       AuszahlungGarantiert = {};      % Nettomonatsauszahlung (garantiert) (tarif, beispiel)(:,szenario)
       AuszahlungFlexibel = {};        % Nettomonatsauszahlung (flexibel) (tarif, beispiel)
       
       % Modellparameter         
       Beleihungswert = 60;
       Horizont = 30;            % Horizont in Jahren fuer Rente
       Aufzinsung = [0.75;       % Hypothetische Rendite fuer Rente
                     1.75; 
                     2.75]      
   
   end

   methods
       
       function obj = rente(zinssatz)
       
          % lade Tarife
          load Versicherung;
          obj.Versicherung = Versicherung;
          clear Versicherung
          
          % zinssatz
          if nargin > 0
              obj.Zinssatz = zinssatz;
          end             
          
          % Rente        
          for anbieter = 1:length(obj.Versicherung)
              for tarif = 1:length(obj.Versicherung(anbieter).RenteGarantiert)
                  
                  % house keeping
                  jahre = (0:obj.Horizont)';
                  aufzinsung = ...
                      [(1 + (obj.Aufzinsung(1)/100)).^jahre,...
                       (1 + (obj.Aufzinsung(2)/100)).^jahre,...
                       (1 + (obj.Aufzinsung(3)/100)).^jahre];
                  
                  % Modellrechnung fuer garantierte Rente                             
                  obj.RenteGarantiert{anbieter,tarif} = ...
                      [obj.Versicherung(anbieter).RenteGarantiert(tarif).*aufzinsung(:,1),...
                      obj.Versicherung(anbieter).RenteGarantiert(tarif).*aufzinsung(:,2),...
                      obj.Versicherung(anbieter).RenteGarantiert(tarif).*aufzinsung(:,3)];                   
                  
                  % Modellrechnung fuer flexible Rente
                  if ~isnan(obj.Versicherung(anbieter).RenteFlexibel(tarif))
                        obj.RenteFlexibel{anbieter,tarif} =...
                            ones(obj.Horizont+1,1) .* obj.Versicherung(anbieter).RenteFlexibel(tarif);
                  else obj.RenteFlexibel{anbieter,tarif} = NaN;
                  end
                  
              end
          end
    
          % Zinsbelastung         
          for anbieter = 1:length(obj.Versicherung)
              for beispiel = 1:length(obj.Versicherung(anbieter).RenteGarantiert)
                  
                obj.Zinsen{anbieter,beispiel} = ones(obj.Horizont+1,1) .* obj.Zinssatz/12 .* obj.Versicherung(anbieter).Einmalzahlung(beispiel);
             
              end
          end
          
          % Nettoauszahlung
          for anbieter = 1:length(obj.Versicherung)
              for beispiel = 1:length(obj.Versicherung(anbieter).RenteGarantiert)
                  
                % Auszahlung Garantiert
                obj.AuszahlungGarantiert{anbieter,beispiel} = ...
                    [obj.RenteGarantiert{anbieter,beispiel}(:,1) - obj.Zinsen{anbieter,beispiel},...
                    obj.RenteGarantiert{anbieter,beispiel}(:,2) - obj.Zinsen{anbieter,beispiel},...
                    obj.RenteGarantiert{anbieter,beispiel}(:,3) - obj.Zinsen{anbieter,beispiel}];
             
                % Auszahlung Flexibel
                obj.AuszahlungFlexibel{anbieter,beispiel} = ...
                    obj.RenteFlexibel{anbieter,beispiel} - obj.Zinsen{anbieter,beispiel};
                                
              end
          end
         
       end % function CONSTRUCTOR
   
   end % methods
end % class
