clc
clear all
stations = readtable('stations.csv');
data2001 = readtable('./csvs_per_year/csvs_per_year/madrid_2001.csv');
data2002 = readtable('./csvs_per_year/csvs_per_year/madrid_2002.csv');
data2003 = readtable('./csvs_per_year/csvs_per_year/madrid_2003.csv');
data2004 = readtable('./csvs_per_year/csvs_per_year/madrid_2004.csv');
data2005 = readtable('./csvs_per_year/csvs_per_year/madrid_2005.csv');
data2006 = readtable('./csvs_per_year/csvs_per_year/madrid_2006.csv');
data2007 = readtable('./csvs_per_year/csvs_per_year/madrid_2007.csv');
data2008 = readtable('./csvs_per_year/csvs_per_year/madrid_2008.csv');
data2009 = readtable('./csvs_per_year/csvs_per_year/madrid_2009.csv');
data2010 = readtable('./csvs_per_year/csvs_per_year/madrid_2010.csv');
data2011 = readtable('./csvs_per_year/csvs_per_year/madrid_2011.csv');
data2012 = readtable('./csvs_per_year/csvs_per_year/madrid_2012.csv');
data2013 = readtable('./csvs_per_year/csvs_per_year/madrid_2013.csv');
data2014 = readtable('./csvs_per_year/csvs_per_year/madrid_2014.csv');
data2015 = readtable('./csvs_per_year/csvs_per_year/madrid_2015.csv');
data2016 = readtable('./csvs_per_year/csvs_per_year/madrid_2016.csv');
data2017 = readtable('./csvs_per_year/csvs_per_year/madrid_2017.csv');

stations2001 = intersect(unique(data2001.station), stations.id, 'rows');
stations2002 = intersect(unique(data2002.station), stations.id, 'rows');
stations2003 = intersect(unique(data2003.station), stations.id, 'rows');
stations2004 = intersect(unique(data2004.station), stations.id, 'rows');
stations2005 = intersect(unique(data2005.station), stations.id, 'rows');
stations2006 = intersect(unique(data2006.station), stations.id, 'rows');
stations2007 = intersect(unique(data2007.station), stations.id, 'rows');
stations2008 = intersect(unique(data2008.station), stations.id, 'rows');
stations2009 = intersect(unique(data2009.station), stations.id, 'rows');
stations2010 = intersect(unique(data2010.station), stations.id, 'rows');
% Dal 2011 in poi tutte le 24 stazioni presenti dentro il csv "stations"
% hanno misurazioni
stations2011 = intersect(unique(data2011.station), stations.id, 'rows');
stations2012 = intersect(unique(data2012.station), stations.id, 'rows');
stations2013 = intersect(unique(data2013.station), stations.id, 'rows');
stations2014 = intersect(unique(data2014.station), stations.id, 'rows');
stations2015 = intersect(unique(data2015.station), stations.id, 'rows');
stations2016 = intersect(unique(data2016.station), stations.id, 'rows');
stations2017 = intersect(unique(data2017.station), stations.id, 'rows');


%% MISURAZIONI INPUT Y
% DSTEM vuole in input le Y come matrice, sulle righe vuole le stazioni,
% sulle colonne i giorni.
% Noi dato che abbiamo misurazioni orarie, usiamo f-HDGM usando come terza
% dimensione le misurazioni orarie ?
% Se invece vogliamo usare HDGM, sulle colonne mettiamo direttamente le
% misurazioni orarie ?

%% COVARIANTE INPUT X
% DSTEM vuole X in matrice a 3 dimensioni con formato AxBxC dove:
% A sono le stazioni
% B sono i regressori
% C sono i giorni






