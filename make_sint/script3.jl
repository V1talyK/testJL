"""Сценарий 3"""
P0 = 80;
Pw_inj  = 350
Pw_prod  = 40

wg1 = 1:8
wg2 = 9:16
wg3 = 17:24
wg4 = 25:32

yst = collect(Iterators.partition(1:length(vd),12))

BHPT[1:32,:] .= 500;
BHPT[33:36,:] .= 1;
#Первый год
for iw=1:8
    BHPT[iw,iw:12] .= 70;
end
#Второй год
BHPT[wg1,yst[2]] .= 50;
BHPT[33,yst[2]] .= Pw_inj - 60;
for iw=9:16
    BHPT[iw,iw+4:24] .= 60;
end
#Третий год
BHPT[wg1,yst[3]] .= Pw_prod+15;
BHPT[wg2,yst[3]] .= Pw_prod+15;

BHPT[33,yst[3]] .= Pw_inj - 50;
BHPT[34,yst[3]] .= Pw_inj - 60;
for iw=17:24
    BHPT[iw,iw+8:36] .= Pw_prod+20;
end
#Четвёртый год
BHPT[wg1,yst[4]] .= Pw_prod+10;
BHPT[wg2,yst[4]] .= Pw_prod+10;
BHPT[wg3,yst[4]] .= Pw_prod+10;
#BHPT[35,yst[4]] .= 150;
BHPT[33,yst[4]] .= Pw_inj - 40;
BHPT[34,yst[4]] .= Pw_inj - 50;

for iw=wg4
    BHPT[iw,iw+12:48] .= Pw_prod+20;
end
#Пятый год
BHPT[wg1,yst[5]] .= Pw_prod+20;
BHPT[wg2,yst[5]] .= Pw_prod+10;
BHPT[wg3,yst[5]] .= Pw_prod+10;
BHPT[wg4,yst[5]] .= Pw_prod+10;

BHPT[33,yst[5]] .= Pw_inj - 30;
BHPT[34,yst[5]] .= Pw_inj - 30;
BHPT[36,yst[5]] .= Pw_inj - 20;

#Шестой год - Седьмой год
for iw in [wg1, wg2, wg3, wg4]
    BHPT[iw,yst[6]] .= 30;
    BHPT[iw,yst[7]] .= 30;
end

for t = 6:7
    BHPT[33,yst[t]] .= Pw_inj-30;
    BHPT[34,yst[t]] .= Pw_inj-20;
    BHPT[36,yst[t]] .= Pw_inj-10;
end

#Восьмой - Десятый год
for ts = yst[8:10]
    for iw in [wg1, wg2, wg3, wg4]
        BHPT[iw,ts] .= 30;
    end
    BHPT[33,ts] .= Pw_inj-60;
    BHPT[34,ts] .= Pw_inj-60;
    BHPT[35,ts] .= Pw_inj-70;
    BHPT[36,ts] .= Pw_inj-60;
end

#Отключения
BHPT[1,49:55] .= 500;
BHPT[3,25:72] .= 500;
BHPT[7,25:84] .= 500;
BHPT[13,37:120] .= 500;
BHPT[25,61:72] .= 500;
BHPT[33,25:26] .= 1;
