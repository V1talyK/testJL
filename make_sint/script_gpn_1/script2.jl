"""Сценарий 1 Газовая шапка"""
P0 = 190;
Pw_inj  = 285
Pw_prod  = 150

wg1 = 1:11      #Добывающие
wg2 = 12:14     #Нагнетательные вод.
wg3 = 15:16     #Нагнетательные газ.

yst = collect(Iterators.partition(1:length(vd),12))

BHPT[wg1,:] .= 500;
BHPT[wg2,:] .= 1;
BHPT[wg3,:] .= 1;
#1-2 год
for iw=wg1
    ts = iw*2:24
    BHPT[iw,ts] .= Pw_prod.+10 .-10 .*mlin(ts);
end
#3-4 год
for iw=1:length(wg1)
    BHPT[wg1[iw],25:48] .= Pw_prod .+wg1[iw] .+ 10. *mlin(25:48);
end
for iw=1:length(wg2)
    BHPT[iw, 24+(iw*4):48] .= Pw_inj;
end
#5-6 год
BHPT[wg1,49:72] .= Pw_prod+15;
BHPT[wg2,49:72] .= Pw_inj+10;
for iw=length(wg3)
    BHPT[iw, 49+(iw*6):72] .= Pw_inj-5;
end

#7-8 год
BHPT[wg1,73:96] .= (Pw_prod.+20 .*mlin(73:96))';
BHPT[wg1[6:end],76:96] .= (Pw_prod.-10 .*mlin(76:96))';
BHPT[wg2,73:96] .= Pw_inj+5;
BHPT[wg3,73:96] .= Pw_inj+10;
#9-10 год
BHPT[wg1,97:120] .= Pw_prod+15;
BHPT[wg2,97:120] .= Pw_inj-5;
BHPT[wg3,97:120] .= Pw_inj+5;
#11-12 год
BHPT[wg1,121:144] .= Pw_prod+10;
BHPT[wg2,121:144] .= Pw_inj-10;
BHPT[wg3,121:144] .= Pw_inj+10;

BHPT[wg1[5],121:144] .= Pw_prod.+10 .-15 .*mlin(121:144);

BHPT[wg1[6],121:132] .= Pw_prod.+15 .-15 .*msofts(121:132);
BHPT[wg1[6],133:144] .= Pw_prod.+0 .-15 .*msofts(133:144);

BHPT[wg1[7],124:132] .= Pw_prod.+10 .-15 .*msofts(124:132);
BHPT[wg1[7],133:144] .= Pw_prod.-5 .-15 .*msofts(133:144);

#13-14 год
BHPT[wg1,145:168] .= Pw_prod+5;
BHPT[wg2,145:168] .= Pw_inj-10;
BHPT[wg3,145:168] .= Pw_inj+5;

BHPT[wg1[5],145:156] .= Pw_prod.+5 .-10 .*msofts(145:156,4,true);
BHPT[wg1[6],145:156] .= Pw_prod.-15 .+15 .*mlin(145:156);
BHPT[wg1[7],145:168] .= Pw_prod.-20 .-15 .*msofts(145:168);
BHPT[wg3[2],145:168] .= Pw_inj+5 .-15 .*msofts(145:168);

#15-16 год
BHPT[wg1,169:192] .= Pw_prod;
BHPT[wg2,169:192] .= Pw_inj;
BHPT[wg3,169:192] .= Pw_inj-5;
#17-20 год
BHPT[wg1,193:240] .= (Pw_prod.+15 .*mlin(193:240))';
BHPT[wg2,193:240] .= Pw_inj;
BHPT[wg3,193:240] .= (Pw_inj.-20 .*mlin(193:240))';
BHPT[wg1[1],193:240] .= Pw_prod.-30 .*msofts(193:240);

#Отключения
BHPT[1,49:55] .= 500;
BHPT[3,25:72] .= 500;
BHPT[7,25:84] .= 500;
BHPT[13,37:120] .= 1;
BHPT[2,120:180] .= 500;
BHPT[13,140:150] .= 1;

BHPT[12,150:170] .= 1;
BHPT[16,200:210] .= 1;

BHPT[4,:] .= 185;
