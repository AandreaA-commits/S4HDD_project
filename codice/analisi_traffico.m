clc
clear all
totaleraw = [];

% ATTENZIONE CHE PER 2023 NON HO OTTOBRE, NOVEMBRE E DICEMBRE

gennaioraw18 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2018.csv');
febbraioraw18 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2018.csv');
marzoraw18 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo2018.csv');
aprileraw18 = readtable('./traffico_madrid_raw/DatosEstacionesAbril2018.csv');
maggioraw18 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2018.csv');
giugnoraw18 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2018.csv');
luglioraw18 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2018.csv');
agostoraw18 = readtable('./traffico_madrid_raw/DatosEstacionesAgosto2018.csv');
settembreraw18 = readtable('./traffico_madrid_raw/DatosEstacionesSeptiembre2018.csv');
ottobreraw18 = readtable('./traffico_madrid_raw/DatosEstacionesOctubre2018.csv');
novembreraw18 = readtable('./traffico_madrid_raw/DatosEstacionesNoviembre2018.csv');
dicembreraw18 = readtable('./traffico_madrid_raw/DatosEstacionesDiciembre2018.csv');

anno2018raw = [];
anno2018raw{1,1}=gennaioraw18;
anno2018raw{1,2}=31;
anno2018raw{2,1}=febbraioraw18;
anno2018raw{2,2}=28;
anno2018raw{3,1}=marzoraw18;
anno2018raw{3,2}=31;
anno2018raw{4,1}=aprileraw18;
anno2018raw{4,2}=30;
anno2018raw{5,1}=maggioraw18;
anno2018raw{5,2}=31;
anno2018raw{6,1}=giugnoraw18;
anno2018raw{6,2}=30;
anno2018raw{7,1}=luglioraw18;
anno2018raw{7,2}=31;
anno2018raw{8,1}=agostoraw18;
anno2018raw{8,2}=31;
anno2018raw{9,1}=settembreraw18;
anno2018raw{9,2}=30;
anno2018raw{10,1}=ottobreraw18;
anno2018raw{10,2}=31;
anno2018raw{11,1}=novembreraw18;
anno2018raw{11,2}=30;
anno2018raw{12,1}=dicembreraw18;
anno2018raw{12,2}=31;


gennaioraw19 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2019.csv');
febbraioraw19 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2019.csv');
marzoraw19 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo2019.csv');
aprileraw19 = readtable('./traffico_madrid_raw/DatosEstacionesAbril2019.csv');
maggioraw19 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2019.csv');
giugnoraw19 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2019.csv');
luglioraw19 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2019.csv');
agostoraw19 = readtable('./traffico_madrid_raw/DatosEstacionesAgosto2019.csv');
settembreraw19 = readtable('./traffico_madrid_raw/DatosEstacionesSeptiembre2019.csv');
ottobreraw19 = readtable('./traffico_madrid_raw/DatosEstacionesOctubre2019.csv');
novembreraw19 = readtable('./traffico_madrid_raw/DatosEstacionesNoviembre2019.csv');
dicembreraw19 = readtable('./traffico_madrid_raw/DatosEstacionesDiciembre2019.csv');

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


gennaioraw20 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2020.csv');
febbraioraw20 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2020.csv');
marzoraw20 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo_2020.csv');
aprileraw20 = readtable('./traffico_madrid_raw/DatosEstacionesAbril2020.csv');
maggioraw20 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2020.csv');
giugnoraw20 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2020.csv');
luglioraw20 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2020.csv');
agostoraw20 = readtable('./traffico_madrid_raw/DatosEstacionesAgosto2020.csv');
settembreraw20 = readtable('./traffico_madrid_raw/DatosEstacionesSeptiembre2020.csv');
ottobreraw20 = readtable('./traffico_madrid_raw/DatosEstacionesOctubre2020.csv');
novembreraw20 = readtable('./traffico_madrid_raw/DatosEstacionesNoviembre2020.csv');
dicembreraw20 = readtable('./traffico_madrid_raw/DatosEstacionesDiciembre2020.csv');

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


gennaioraw21 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2021.csv');
febbraioraw21 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2021.csv');
marzoraw21 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo2021.csv');
aprileraw21 = readtable('./traffico_madrid_raw/DatosEstacionesAbril2021.csv');
maggioraw21 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2021.csv');
giugnoraw21 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2021.csv');
luglioraw21 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2021.csv');
agostoraw21 = readtable('./traffico_madrid_raw/DatosEstacionesAGOSTO_2021.csv');
settembreraw21 = readtable('./traffico_madrid_raw/DatosEstacionesSEPTIEMBRE_2021.csv');
ottobreraw21 = readtable('./traffico_madrid_raw/DatosEstacionesOCTUBRE_2021.csv');
novembreraw21 = readtable('./traffico_madrid_raw/DatosEstacionesNOVIEMBRE_2021.csv');
dicembreraw21 = readtable('./traffico_madrid_raw/DatosEstacionesDICIEMBRE_2021.csv');

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


gennaioraw22 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2022.csv');
febbraioraw22 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2022.csv');
marzoraw22 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo2022.csv');
aprileraw22 = readtable('./traffico_madrid_raw/DatosEstacionesAbril2022.csv');
maggioraw22 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2022.csv');
giugnoraw22 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2022.csv');
luglioraw22 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2022.csv');
agostoraw22 = readtable('./traffico_madrid_raw/DatosEstacionesAgosto2022.csv');
settembreraw22 = readtable('./traffico_madrid_raw/DatosEstacionesSeptiembre2022.csv');
ottobreraw22 = readtable('./traffico_madrid_raw/DatosEstacionesOctubre2022.csv');
novembreraw22 = readtable('./traffico_madrid_raw/DatosEstacionesNoviembre2022.csv');
dicembreraw22 = readtable('./traffico_madrid_raw/DatosEstacionesDiciembre2022.csv');

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


gennaioraw23 = readtable('./traffico_madrid_raw/DatosEstacionesEnero2023.csv');
febbraioraw23 = readtable('./traffico_madrid_raw/DatosEstacionesFebrero2023.csv');
marzoraw23 = readtable('./traffico_madrid_raw/DatosEstacionesMarzo2023.csv');
aprileraw23 = readtable('./traffico_madrid_raw/DatosEstacionesAbril023.csv');
maggioraw23 = readtable('./traffico_madrid_raw/DatosEstacionesMayo2023.csv');
giugnoraw23 = readtable('./traffico_madrid_raw/DatosEstacionesJunio2023.csv');
luglioraw23 = readtable('./traffico_madrid_raw/DatosEstacionesJulio2023.csv');
agostoraw23 = readtable('./traffico_madrid_raw/DatosEstacionesagosto2023.csv');
settembreraw23 = readtable('./traffico_madrid_raw/DatosEstacionesSeptiembre2023.csv');

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


totaleraw{1,1} = anno2018raw;
totaleraw{2,1} = anno2019raw;
totaleraw{3,1} = anno2020raw;
totaleraw{4,1} = anno2021raw;
totaleraw{5,1} = anno2022raw;
totaleraw{6,1} = anno2023raw;

%%

totale = [];
for idxanno = 1:6
    anno = [];
    for idxmese = 1:height(totaleraw{idxanno,1})
        data = totaleraw{idxanno,1}{idxmese,1};
        giorni = totaleraw{idxanno,1}{idxmese,2}    
        from = 1;
        to = 4;
        mese = [];
        for day = 1:giorni
            ap = [];
            for station = 1:59
                ap18 = data(from : to, 4:15);
                fin18 = [ap18(1,:)+ap18(3,:); ap18(2,:)+ap18(4,:)];
                ff = reshape(table2array(fin18),1,[]);
                ap = [ap;ff];
                from = to + 1;
                to = from + 3;
            end
            mese = [mese ap];
        end
        anno = [anno mese];
    end
    totale = [totale anno];
end
%%
writematrix(totale,'dati_traffico_ready.csv');


