clc
clear all

dati_meteo_2019raw = [];

% ATTENZIONE CHE PER 2023 NON HO OTTOBRE, NOVEMBRE E DICEMBRE

gennaioraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\a_ene_mo19.csv');
febbraioraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\b_feb_mo19.csv');
marzoraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\c_mar_mo19.csv');
aprileraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\d_abr_mo19.csv');
maggioraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\e_may_mo19.csv');
giugnoraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\f_jun_mo19.csv');
luglioraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\g_jul_mo19.csv');
agostoraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\h_ago_mo19.csv');
settembreraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\i_sep_mo19.csv');
ottobreraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\l_oct_mo19.csv');
novembreraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\m_nov_mo19.csv');
dicembreraw19 = readtable('C:\Users\arici\Desktop\scraper_gabs\Anio201912\n_dic_mo19.csv');

dati_staz = readtable('C:\Users\arici\Desktop\scriping dati\csv_dataset\stations.csv');

anno2019raw = [];
anno2019raw{1,1}=gennaioraw19;
anno2019raw{1,2}=31;
anno2019raw{2,1}=febbraioraw19;
anno2019raw{2,2}=28;
anno2019raw{3,1}=marzoraw19;
anno2019raw{3,2}=31;
anno2019raw{4,1}=aprileraw19;
anno2019raw{4,2}=30;
anno2019raw{5,1}=maggioraw19;
anno2019raw{5,2}=31;
anno2019raw{6,1}=giugnoraw19;
anno2019raw{6,2}=30;
anno2019raw{7,1}=luglioraw19;
anno2019raw{7,2}=31;
anno2019raw{8,1}=agostoraw19;
anno2019raw{8,2}=31;
anno2019raw{9,1}=settembreraw19;
anno2019raw{9,2}=30;
anno2019raw{10,1}=ottobreraw19;
anno2019raw{10,2}=31;
anno2019raw{11,1}=novembreraw19;
anno2019raw{11,2}=30;
anno2019raw{12,1}=dicembreraw19;
anno2019raw{12,2}=31;

dati_meteo_2019raw{1,1} = anno2019raw;

%abbiamo una table per ogni mese
totale = [];

giorni_del_mese = [31 28 31 30 31 30 31 31 30 31 30 31];

mesi_disponibili = [12,12,12,12,9];

for anno = 1:1
    dati_anno = dati_meteo_2019raw{anno,1};
    stazioni = unique(dati_anno{1,1}.ESTACION);
    covariate = unique(dati_anno{1,1}.MAGNITUD);        
    for mese = 1:size(dati_anno,1)
        dati_mese = dati_anno{mese,1};
        %stazioni = unique(dati_mese.ESTACION);
        %covariate = unique(dati_mese.MAGNITUD);
        i = 1; % index covariata
        for meteo_covariata = covariate'
            covariata_filtered = dati_mese(find(dati_mese.MAGNITUD == meteo_covariata), :);
            lista_stazioni_singola_covariata = [];
            for station = stazioni'
                station_filtered = covariata_filtered(find(covariata_filtered.ESTACION == station),:); 
                stazione = [];     
                offsetsx = 0;  
                offsetdx = 0;  
                for row = 1:size(station_filtered, 1)
                    if isempty(stazione)
                        stazione = [station stazione];
                        f = find(mod(dati_staz{:,1},100) == station);
                        % id elevazione latitudine longitudine osservazioni
                        stazione = [stazione [dati_staz{f,6}, dati_staz{f,5}, dati_staz{f,4}]];
                    end
                    aggiunta = station_filtered(row, 9:end);
                    % nel caso manchi un giorno in cima alle righe, aggiungo una riga di NaN
                    if station_filtered(row, :).DIA+offsetsx > row+offsetdx
                        disp("giorno mancante dx");
                        stazione = [stazione NaN(1,24)];
                        offsetdx = offsetdx + 1;
                    % nel caso manchi un giorno in mezzo alle righe
                    elseif station_filtered(row, :).DIA+offsetsx < row+offsetdx
                        disp("giorno mancante sx");
                        stazione = [stazione NaN(1,24)];
                        offsetsx = offsetsx + 1;
                    end
                    for j = 1:2:size(aggiunta, 2)
                        if cell2mat(aggiunta{:, j+1}) == 'N'
                            aggiunta{:,j} = NaN(1,1);
                        end
                    end                    
                    stazione = [stazione aggiunta{:, 1:2:end}]; 
                    % nel caso manchino giorni in fondo alle righe
                    if size(station_filtered,1) < (giorni_del_mese(:,mese) - (offsetsx+offsetdx)) & size(station_filtered, 1) ~= 0 & row == size(station_filtered,1)
                        for mm = 1:(giorni_del_mese(:,mese) - size(station_filtered,1) - (offsetsx+offsetdx))
                            stazione  = [stazione NaN(1,24)];
                        end
                        %stazione = [stazione NaN(1, 24*(giorni_del_mese(:,mese) - size(station_filtered,1) - (offsetsx+offsetdx)))];
                    end    
                end                              
                lista_stazioni_singola_covariata = [lista_stazioni_singola_covariata; stazione];
            end 
            dati_meteo_2019{i,1} = meteo_covariata;            
            dati_meteo_2019{i,mese+1} = array2table(lista_stazioni_singola_covariata);              
            i = i+1;    
        end
    end
end


tabella_poll = table();

tabella_poll.CODICE_POLLUTANT = [1 6 7 8 9 10 12 14 20 30 35 37 38 39 42 43 44]';
tabella_poll.NOME = {'SO2' 'CO' 'NO' 'NO2' 'PM2.5' 'PM10' 'NOX' 'O3' 'TOL' 'BEN' 'EBE' 'MXY' 'PXY' 'OXY' 'TCH' 'CH4' 'NMHC'}';


%sostituzione del codice della stazione con lalatitudine e longitudine

for i=1:size(dati_meteo_2019,1)
    indice = find(tabella_poll.CODICE_POLLUTANT == dati_meteo_2019{i,1});
    dati_meteo_2019{i,1} = tabella_poll.NOME(indice);
end

save("dati_pollutants_2019", "dati_meteo_2019");


