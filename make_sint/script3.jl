"""Сценарий 3"""
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
BHPT[33,yst[2]] .= 120;
for iw=9:16
    BHPT[iw,iw+4:24] .= 60;
end
#Третий год
BHPT[wg1,yst[3]] .= 40;
BHPT[wg2,yst[3]] .= 40;

BHPT[33,yst[3]] .= 130;
BHPT[34,yst[3]] .= 120;
for iw=17:24
    BHPT[iw,iw+8:36] .= 50;
end
#Четвёртый год
BHPT[wg1,yst[4]] .= 30;
BHPT[wg2,yst[4]] .= 30;
BHPT[wg3,yst[4]] .= 30;
#BHPT[35,yst[4]] .= 150;
BHPT[33,yst[4]] .= 140;
BHPT[34,yst[4]] .= 130;

for iw=wg4
    BHPT[iw,iw+12:48] .= 40;
end
#Пятый год
BHPT[wg1,yst[5]] .= 40;
BHPT[wg2,yst[5]] .= 30;
BHPT[wg3,yst[5]] .= 30;
BHPT[wg4,yst[5]] .= 30;

BHPT[33,yst[5]] .= 150;
BHPT[34,yst[5]] .= 150;
BHPT[36,yst[5]] .= 160;

#Шестой год - Седьмой год
for iw in [wg1, wg2, wg3, wg4]
    BHPT[iw,yst[6]] .= 30;
    BHPT[iw,yst[7]] .= 30;
end

for t = 6:7
    BHPT[33,yst[t]] .= 120;
    BHPT[34,yst[t]] .= 110;
    BHPT[36,yst[t]] .= 100;
end

#Восьмой - Десятый год
for ts = yst[8:10]
    for iw in [wg1, wg2, wg3, wg4]
        BHPT[iw,ts] .= 30;
    end
    BHPT[33,ts] .= 90;
    BHPT[34,ts] .= 90;
    BHPT[35,ts] .= 80;
    BHPT[36,ts] .= 90;
end

#Отключения
BHPT[1,49:55] .= 500;
BHPT[25,61:72] .= 500;
BHPT[7,25:84] .= 500;
BHPT[33,25:26] .= 1;
