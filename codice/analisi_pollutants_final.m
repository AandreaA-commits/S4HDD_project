clc
clear all

dati_pollutants_2019raw = [];

% ATTENZIONE CHE PER 2023 NON HO OTTOBRE, NOVEMBRE E DICEMBRE

gennaioraw19 = readtable('./pollutants_madrid_raw/a_ene_mo19.csv');
febbraioraw19 = readtable('./pollutants_madrid_raw/b_feb_mo19.csv');
marzoraw19 = readtable('./pollutants_madrid_raw/c_mar_mo19.csv');
aprileraw19 = readtable('./pollutants_madrid_raw/d_abr_mo19.csv');
maggioraw19 = readtable('./pollutants_madrid_raw/e_may_mo19.csv');
giugnoraw19 = readtable('./pollutants_madrid_raw/f_jun_mo19.csv');
luglioraw19 = readtable('./pollutants_madrid_raw/g_jul_mo19.csv');
agostoraw19 = readtable('./pollutants_madrid_raw/h_ago_mo19.csv');
settembreraw19 = readtable('./pollutants_madrid_raw/i_sep_mo19.csv');
ottobreraw19 = readtable('./pollutants_madrid_raw/l_oct_mo19.csv');
novembreraw19 = readtable('./pollutants_madrid_raw/m_nov_mo19.csv');
dicembreraw19 = readtable('./pollutants_madrid_raw/n_dic_mo19.csv');

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

%{
gennaioraw20 = readtable('./pollutants_madrid_raw/ene_meteo20.csv');
febbraioraw20 = readtable('./pollutants_madrid_raw/feb_meteo20.csv');
marzoraw20 = readtable('./pollutants_madrid_raw/mar_meteo20.csv');
aprileraw20 = readtable('./pollutants_madrid_raw/abr_meteo20.csv');
maggioraw20 = readtable('./pollutants_madrid_raw/may_meteo20.csv');
giugnoraw20 = readtable('./pollutants_madrid_raw/jun_meteo20.csv');
luglioraw20 = readtable('./pollutants_madrid_raw/jul_meteo20.csv');
agostoraw20 = readtable('./pollutants_madrid_raw/ago_meteo20.csv');
settembreraw20 = readtable('./pollutants_madrid_raw/sep_meteo20.csv');
ottobreraw20 = readtable('./pollutants_madrid_raw/oct_meteo20.csv');
novembreraw20 = readtable('./pollutants_madrid_raw/nov_meteo20.csv');
dicembreraw20 = readtable('./pollutants_madrid_raw/dic_meteo20.csv');

anno2020raw = [];
anno2020raw{1,1}=gennaioraw20;
anno2020raw{1,2}=31;
anno2020raw{2,1}=febbraioraw20;
anno2020raw{2,2}=28;
anno2020raw{3,1}=marzoraw20;
anno2020raw{3,2}=31;
anno2020raw{4,1}=aprileraw20;
anno2020raw{4,2}=30;
anno2020raw{5,1}=maggioraw20;
anno2020raw{5,2}=31;
anno2020raw{6,1}=giugnoraw20;
anno2020raw{6,2}=30;
anno2020raw{7,1}=luglioraw20;
anno2020raw{7,2}=31;
anno2020raw{8,1}=agostoraw20;
anno2020raw{8,2}=31;
anno2020raw{9,1}=settembreraw20;
anno2020raw{9,2}=30;
anno2020raw{10,1}=ottobreraw20;
anno2020raw{10,2}=31;
anno2020raw{11,1}=novembreraw20;
anno2020raw{11,2}=30;
anno2020raw{12,1}=dicembreraw20;
anno2020raw{12,2}=31;


gennaioraw21 = readtable('./pollutants_madrid_raw/ene_meteo21.csv');
febbraioraw21 = readtable('./pollutants_madrid_raw/feb_meteo21.csv');
marzoraw21 = readtable('./pollutants_madrid_raw/mar_meteo21.csv');
aprileraw21 = readtable('./pollutants_madrid_raw/abr_meteo21.csv');
maggioraw21 = readtable('./pollutants_madrid_raw/may_meteo21.csv');
giugnoraw21 = readtable('./pollutants_madrid_raw/jun_meteo21.csv');
luglioraw21 = readtable('./pollutants_madrid_raw/jul_meteo21.csv');
agostoraw21 = readtable('./pollutants_madrid_raw/ago_meteo21.csv');
settembreraw21 = readtable('./pollutants_madrid_raw/sep_meteo21.csv');
ottobreraw21 = readtable('./pollutants_madrid_raw/oct_meteo21.csv');
novembreraw21 = readtable('./pollutants_madrid_raw/nov_meteo21.csv');
dicembreraw21 = readtable('./pollutants_madrid_raw/dic_meteo21.csv');

anno2021raw = [];
anno2021raw{1,1}=gennaioraw21;
anno2021raw{1,2}=31;
anno2021raw{2,1}=febbraioraw21;
anno2021raw{2,2}=28;
anno2021raw{3,1}=marzoraw21;
anno2021raw{3,2}=31;
anno2021raw{4,1}=aprileraw21;
anno2021raw{4,2}=30;
anno2021raw{5,1}=maggioraw21;
anno2021raw{5,2}=31;
anno2021raw{6,1}=giugnoraw21;
anno2021raw{6,2}=30;
anno2021raw{7,1}=luglioraw21;
anno2021raw{7,2}=31;
anno2021raw{8,1}=agostoraw21;
anno2021raw{8,2}=31;
anno2021raw{9,1}=settembreraw21;
anno2021raw{9,2}=30;
anno2021raw{10,1}=ottobreraw21;
anno2021raw{10,2}=31;
anno2021raw{11,1}=novembreraw21;
anno2021raw{11,2}=30;
anno2021raw{12,1}=dicembreraw21;
anno2021raw{12,2}=31;


gennaioraw22 = readtable('./pollutants_madrid_raw/ene_meteo22.csv');
febbraioraw22 = readtable('./pollutants_madrid_raw/feb_meteo22.csv');
marzoraw22 = readtable('./pollutants_madrid_raw/mar_meteo22.csv');
aprileraw22 = readtable('./pollutants_madrid_raw/abr_meteo22.csv');
maggioraw22 = readtable('./pollutants_madrid_raw/may_meteo22.csv');
giugnoraw22 = readtable('./pollutants_madrid_raw/jun_meteo22.csv');
luglioraw22 = readtable('./pollutants_madrid_raw/jul_meteo22.csv');
agostoraw22 = readtable('./pollutants_madrid_raw/ago_meteo22.csv');
settembreraw22 = readtable('./pollutants_madrid_raw/sep_meteo22.csv');
ottobreraw22 = readtable('./pollutants_madrid_raw/oct_meteo22.csv');
novembreraw22 = readtable('./pollutants_madrid_raw/nov_meteo22.csv');
dicembreraw22 = readtable('./pollutants_madrid_raw/dic_meteo22.csv');

anno2022raw = [];
anno2022raw{1,1}=gennaioraw22;
anno2022raw{1,2}=31;
anno2022raw{2,1}=febbraioraw22;
anno2022raw{2,2}=28;
anno2022raw{3,1}=marzoraw22;
anno2022raw{3,2}=31;
anno2022raw{4,1}=aprileraw22;
anno2022raw{4,2}=30;
anno2022raw{5,1}=maggioraw22;
anno2022raw{5,2}=31;
anno2022raw{6,1}=giugnoraw22;
anno2022raw{6,2}=30;
anno2022raw{7,1}=luglioraw22;
anno2022raw{7,2}=31;
anno2022raw{8,1}=agostoraw22;
anno2022raw{8,2}=31;
anno2022raw{9,1}=settembreraw22;
anno2022raw{9,2}=30;
anno2022raw{10,1}=ottobreraw22;
anno2022raw{10,2}=31;
anno2022raw{11,1}=novembreraw22;
anno2022raw{11,2}=30;
anno2022raw{12,1}=dicembreraw22;
anno2022raw{12,2}=31;


gennaioraw23 = readtable('./pollutants_madrid_raw/ene_meteo23.csv');
febbraioraw23 = readtable('./pollutants_madrid_raw/feb_meteo23.csv');
marzoraw23 = readtable('./pollutants_madrid_raw/mar_meteo23.csv');
aprileraw23 = readtable('./pollutants_madrid_raw/abr_meteo23.csv');
maggioraw23 = readtable('./pollutants_madrid_raw/may_meteo23.csv');
giugnoraw23 = readtable('./pollutants_madrid_raw/jun_meteo23.csv');
luglioraw23 = readtable('./pollutants_madrid_raw/jul_meteo23.csv');
agostoraw23 = readtable('./pollutants_madrid_raw/ago_meteo23.csv');
settembreraw23 = readtable('./pollutants_madrid_raw/sep_meteo23.csv');

anno2023raw = [];
anno2023raw{1,1}=gennaioraw23;
anno2023raw{1,2}=31;
anno2023raw{2,1}=febbraioraw23;
anno2023raw{2,2}=28;
anno2023raw{3,1}=marzoraw23;
anno2023raw{3,2}=31;
anno2023raw{4,1}=aprileraw23;
anno2023raw{4,2}=30;
anno2023raw{5,1}=maggioraw23;
anno2023raw{5,2}=31;
anno2023raw{6,1}=giugnoraw23;
anno2023raw{6,2}=30;
anno2023raw{7,1}=luglioraw23;
anno2023raw{7,2}=31;
anno2023raw{8,1}=agostoraw23;
anno2023raw{8,2}=31;
anno2023raw{9,1}=settembreraw23;
anno2023raw{9,2}=30;
%}

dati_pollutants_2019raw{1,1} = anno2019raw;
%dati_pollutants_2019raw{2,1} = anno2020raw;
%dati_pollutants_2019raw{3,1} = anno2021raw;
%dati_pollutants_2019raw{4,1} = anno2022raw;
%dati_pollutants_2019raw{5,1} = anno2023raw;

%abbiamo una table per ogni mese
dati_pollutants_2019 = [];

giorni_del_mese = [31 28 31 30 31 30 31 31 30 31 30 31];

mesi_disponibili = [12,12,12,12,9];

dati_staz = readtable('./stazioni_raw/dati_stazioni_pollutants.csv');


for anno = 1:1
    dati_anno = dati_pollutants_2019raw{anno,1};
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
                        f = find(dati_staz{:,1} == (28079000 + station));
                        stazione = [stazione [dati_staz{f,6}, dati_staz{f,5},dati_staz{f,4}]];
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
            dati_pollutants_2019{i,1} = meteo_covariata;            
            dati_pollutants_2019{i,mese+1} = array2table(lista_stazioni_singola_covariata);              
            i = i+1;    
        end
    end
end


for i = 1:12
    unique(dati_pollutants_2019raw{1,1}{i,1}.MAGNITUD)
end

tabella_poll = table();

tabella_poll.CODICE_POLLUTANT = [1 6 7 8 9 10 12 14 20 30 35 37 38 39 42 43 44]';
tabella_poll.NOME = {'SO2' 'CO' 'NO' 'NO2' 'PM2.5' 'PM10' 'NOX' 'O3' 'TOL' 'BEN' 'EBE' 'MXY' 'PXY' 'OXY' 'TCH' 'CH4' 'NMHC'}';


%sostituzione del codice della stazione con lalatitudine e longitudine

for i=1:size(dati_pollutants_2019,1)
    indice = find(tabella_poll.CODICE_POLLUTANT == dati_pollutants_2019{i,1});
    dati_pollutants_2019{i,1} = tabella_poll.NOME(indice);
end

save ("dati_pollutants_2019", "dati_pollutants_2019");
