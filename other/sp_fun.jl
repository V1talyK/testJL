function Solve_k(aak, k_giv, xkk, ykk)
    ### Загрузить входные данные
    itermax1 = 600; iterkmax = 81; i_int =1; Lambda1 = 0.0251;
    
    Nkk1 = 9;  Nkk2 = 9;#   Nkk =length(xkk); 
    Nkk = Nkk1*Nkk2;
    xkk = zeros(Nkk,1); ykk = xkk;
    
    for i=1:Nkk1
        for j=1:Nkk2
            ij = (i-1)*Nkk2+j;
            xkk[ij] = (i-1)*2;
            ykk[ij] = (j-1)*2;
        end # Nkk
    end # Nkk
    
    # k_giv = 10.0*0.00864*ones(Nkk,1);  #  k_giv = load('k_giv.txt');
     k_giv = (2 + xkk + ykk);
     scatter3(xkk',ykk',k_giv/0.00864);
     aak_tr = zeros(Nkk,1);
     aak_te = ones(Nkk,1);
    
    ii[1] = 5; ii[2] = 15; ii[3] = 55; ii[4] = 42; ii[4] = 25; ii[5] = 70;
    for i=1:length(ii)
        aak_tr[ii[i]] = 1; aak_te[ii[i]] = 0;
    end # Nkk
    
    #  Визуализация точек с номерами
    #graf1m([xkk',ykk']);
    graf1m([xkk,ykk]);
    ### Задать начальное приближение для решения обратной задачи и базисных функций
    
     k0kk = 10.0*0.00864*ones(Nkk,1); #fff = 1*55*0.00864*rand(1,length(k0kk)); k0kk = k0kk + fff';
    # k0kk = k_giv;
     scatter3(xkk',ykk',k0kk/0.00864);
      # trisurf(delaunay(xkk,ykk),xkk,ykk,k0kk/0.00864);
    
    Fim =  ones(Nkk,1); # Fim  = spdiags(ones(Nkk,1),0,Nkk,Nkk);  #   spy(Fim)
    
    A_k0kk = zeros(Nkk,iterkmax);  J1iterk = zeros(Nkk,1); J2iterk = zeros(Nkk,1);
        
    xkkm = (xkk - min(xkk))/(max(xkk) - min(xkk));
    ykkm = (ykk - min(ykk))/(max(ykk) - min(ykk));
    xkk1 = xkkm';
    ykk1 = ykkm';
    ######################## Начало цикла выделения кластеров. Идентификация Fim ####################
    
    iterk =0;
    while (iterk < iterkmax)
       iterk = iterk + 1;
    
    ######## Решение обратной задачи с заданной матрицей Fim и расчет целевого функционала на обучающей выборке #######
            
        [k0e,JJ1,Dww] = Calc_kmod_i_1(itermax1,k0kk,k_giv,aak_tr,Fim,i_int,Lambda1);
        
        k0kk =  k0e;   scatter3(xkk1,ykk1,k0kk/0.00864);
        J1iterk[iterk] = JJ1;
        A_k0kk[:,iterk] = k0kk;
        ####################  Расчет целевого функционала на тестовой выборке
        
        JJ2 = Obj_function(k0kk,k_giv,aak_te);
        J2iterk[iterk] = JJ2;
        
        ####################  Расчет новых значений Fim
        Dwwm = (Dww - min(Dww))/(max(Dww) - min(Dww));  #    Dww = Fim*(((Dww'*Fim)./(ones(1,Nkk)*Fim))');  #scatter3(xkk1,ykk1,Dwwm)
        
        #  выделение кластеров
        X2 = [1.0*xkk1',1.0*ykk1',Dwwm]; 
        Z = linkage(X2,"ward");
        #  формирование Fim через выбранные кластеры (функцию cluster)
        ckk = cluster(Z,"Maxclust",iterk+1);
        Nclust = max(ckk);
        Fim = sparse([],[],[],Nkk,Nclust);

        for kk =1:Nkk
            kklust = ckk[kk];
            Fim[kk,kklust] = 1;
        end
    
     #  Fim = ones(Nkk,1);  #   Fim = spdiags(ones(Nkk,1),0,Nkk,Nkk);   # spy(Fim)
    
    end    #while(Iterk)
    
    scatter3(xkk',ykk',k0kk/0.00864);
    
    ####################  Конец выделения кластеров    ###
    # Jk = zeros(Nkk,1);
    #  ### ### ###   ###             Сохранить результаты
    #    ###  #
    #   save('xkkt.txt','xkk','-ascii');  save('ykkt.txt','ykk','-ascii'); save('A_k0kkt.txt','A_k0kk','-ascii')  ;
    #       ee =  load('A_k0kkt.txt');
    #       k0kk = ee(:,5);
    #      scatter3(xkk',ykk',k0kk/0.00864);
    #    trisurf(delaunay(xkk,ykk),xkk,ykk,k0kk/0.00864);
    
    
    #       save('J1iterk.txt','J1iterk','-ascii'); save('J2iterk.txt','J2iterk','-ascii'); save('Jk.txt','Jk','-ascii');
    #  # визуализация целевого функционала
    #  hgh1 = J1iterk/J1iterk(1);  hgh2 = J2iterk/J2iterk(1);
    #   plot(-log(hgh1),-log(hgh2))
    #     hold on
    #  plot(-log(hgh1),-log(hgh2),'*')
    #   hold off
    # #
    #    plot(J1iterk)
    #    hold on
    #  plot(J1iterk,'*')
    #   hold off
    #   #
    #      plot(J2iterk)
    #    hold on
    #  plot(J2iterk,'*')
    #   hold off
    # #
    #      plot(Jk,'*')
    #    hold on
    #       plot(Jk)
    
    
    return k0,alfaW,J
end
    
function Obj_function(Fun,Fun_giv,aa)
    dF = Fun-Fun_giv;
    JJ =  sum(sum(dF.*(aa.*dF)));
    return JJ
end
    
function Clust(I,Z)

#  определение элементов в кластере i

    NN = 1+ length(Z[:,1]); ZZ =  Z[I,:]; aa = [ZZ[1],ZZ[2]]; cc = [];

    iter =0;
    while (iter<1)

    aa1 = [];
        for ii =1:length(aa)
            if (aa[ii]>NN)
                I =  aa(ii)-NN;
                aa1 = [aa1,Z[I,1],Z[I,2]];
            else
                cc = [cc,aa[ii]];
            end
        end

    aa = aa1;

        if (length(aa1)<1)
            iter =2;
        end
    end

    c1 =  ones(NN,1);
    #     cc1 = Clust_old(I,Z);
    cc1 = cc;
    for kk1 =1:length(cc1)
        ii = cc1[kk1];
        c1[ii] = 2;
    end
    return c1
end
    
    
    
function Clust_old(I,Z)
#  определение элементов в кластере i
# I = 40;
    NN = 1+ length(Z[:,1]); ZZ =  Z[I,:]; aa = [ZZ[1],ZZ[2]]; cc = [];
    # Nclust = NN-I;
    iter =0;
    while (iter<1)

    aa1 = [];
    for ii =1:length(aa)
        if (aa(ii)>NN)
            I =  aa(ii)-NN;
            aa1 = [aa1,Z(I,1),Z(I,2)];
        else
            cc = [cc,aa(ii)];
        end
#
    end

    aa = aa1;

    if(length(aa1)<1)
        iter =2;
    end

    end

    return cc
end
    
    
    
    
# function graf1m(X)
#     figure()
# #X = [-1.5 3.2; 1.8 3.3; -3.7 1.5; -1.5 1.3; ...
# #     0.8 1.2; 3.3 1.5; -4.0 -1.0;-2.3 -0.7; ...
# #     0 -0.5; 2.0 -1.5; 3.7 -0.8; -3.5 -2.9; ...
# #    -0.9 -3.9; 2.0 -3.5; 3.5 -2.25];

# # voronoi(X(:,1),X(:,2))
# #  scatter(X(:,1),X(:,2),10)
# scatter(X(:,1),X(:,2),20,'filled')
# box on;
# hold on
# title('Fig.2b');
# hold on


# # Assign labels to the points.
# nump = size(X,1);
# plabels = arrayfun(@(n) {sprintf('# d', n)}, (1:nump)');
# # hold on
# # Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', ...
# #       'bold', 'HorizontalAlignment','left', ...
# #       'BackgroundColor', 'none');

# Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', 'HorizontalAlignment','left', 'BackgroundColor', 'none');

# end
    
    
    
function graf1m1(X)
    
    figure()
    #X = [-1.5 3.2; 1.8 3.3; -3.7 1.5; -1.5 1.3; ...
     #     0.8 1.2; 3.3 1.5; -4.0 -1.0;-2.3 -0.7; ...
     #     0 -0.5; 2.0 -1.5; 3.7 -0.8; -3.5 -2.9; ...
     #    -0.9 -3.9; 2.0 -3.5; 3.5 -2.25];
    
    # voronoi(X(:,1),X(:,2))
    
     scatter(X(:,1),X(:,2),10)
    
    
    # Assign labels to the points.
    nump = size(X,1);
    plabels = arrayfun(@(n) {sprintf('#d', n)}, (1:nump)');

    Hpl = text(X(:,1), X(:,2), plabels, 'FontWeight', 'bold', 'HorizontalAlignment','left', 'BackgroundColor', 'none');
 
end
    
    
    
function Calc_kmod_i_1(itermax,k0kk,k_giv,aa,Fim,i_int,Lambda1)
    
    iter =0;  Jold =1.e+30; J = Jold; III =0;   Nkk =length(k_giv);
    while (iter<itermax)
        iter = iter + 1;
        J = Obj_function(k0kk,k_giv,aa);
    ############ Расчет градиента целевого функционала + Коррекция управляющих параметров (U=a)
    
        Dww2 = 2*aa.*(k0kk-k_giv);
    
      #  Dww2 = Grad_kk(ia,ja,iaw,jaw,iab,jab,alfaW,alfwL,alfb,k0kk,k_giv,aak,alf1,alf2,Sbon,typebon, Ptim, Pw, DQ, PSItim,L,Lb,ii,jj,ibb,iww);
    
        Dww = Fim*(((Dww2'*Fim)./(ones(1,Nkk)*Fim))');
    
        k0kk = k0kk.*(ones(length(k0kk),1) - (Lambda1/(1.e-30 + max(abs(Dww))))*Dww);
    
        dJ = J-Jold;
            if dJ>0
                Lambda1 =  0.85*Lambda1;
            end
    
             if dJ >=0
                 III = III + 1;
        #           Lambda1 =  1.5*Lambda1;
             end
    
             if abs(dJ)/Jold < 0.01
             #      Lambda1 =  1.015*Lambda1;
             end
             Jold = J;
    
      #    trisurf(delaunay(xc,yc),xc,yc,k0/0.00864);
      #  trisurf(delaunay(xkk,ykk),xkk,ykk,k0kk/0.00864);
       #    trisurf(delaunay(xkk,ykk),xkk,ykk,Dww/0.00864);
      #    trisurf(delaunay(xkk,ykk),xkk,ykk,Dww2/0.00864);
    
    ####################  Проверка критерия сходимости  if( |  J < 0.01)
        if III>81 || J < 1.e-4
            break
        end
             #
    
      iter
      III
      J
     # lam
    
     end  #while(Iter)
       #   trisurf(delaunay(xc,yc),xc,yc,k0/0.00864);
      #  trisurf(delaunay(xkk,ykk),xkk,ykk,k0kk/0.00864);
    
    k01 = k0kk; J1 = J; Dww1 = Dww2;
    return k01,J1,Dww1 
end