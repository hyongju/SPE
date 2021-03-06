function [voronoi_rg,vornb,vornb2] = p2Voronoi(pos,bnd_pnts)
% clear all;close all;
tol =1e-15;

%% =======================================================
% Order-2 Voronoi Diagram with set of points in 2D/3D polygon
% ========================================================
% version 1.01
% by Hyongju Park
%---------------------------------------------------------
% inputs: bnd_pnts      boundary points                m x 2
%         pos           points inside the boundary     n x 2
%---------------------------------------------------------
% outputs: voronoi_rg   order-2 Voronoi regions        ? x n  
%          vornb        Voronoi neighbors              1 x n
%          vornb2       Voronoi neighbors              1 x ? (repeating
%          neighbors removed)
% =========================================================================
% This functions works for d = 2, 3
% -------------------------------------------------------------------------
% This function requires:
%       vert2lcon.m (Matt Jacobson / Michael Keder)
%       pbisec.m (by me)
%       con2vert.m (Michael Keder)
%       inhull.m (John D'Errico)
% -------------------------------------------------------------------------
% Written by Hyongju Park, hyongju@gmail.com / park334@illinois.edu
% /hjcpark@umich.edu
% Change logs:
% 11 Nov 2016: fixed overlapping cells issues
% 11 Aug 2015: skip error messages (version 1.01) 
% 5  May 2015: initial release (version 1.0)
% =========================================================================
% Known issues:
% Input points must satisfy assumptions such as non co-circularity and
% general position assumptions
% -------------------------------------------------------------------------
%%
[vornb,~,A_org,b_org] = pVoronoi(pos,bnd_pnts);      % 
%%
% h0 = figure('position',[0 0 700 700],'Color',[1 1 1]);
% k=0;
% % for i = 1:size(voronoi_rg{t},1)*size(voronoi_rg{t},2)
% %     col(i,:)= rand(1,3);
% % %     col(i,:) = [i/(size(voronoi_rg{t},1)*size(voronoi_rg{t},2)) 1 1];
% % end
% % col = distinguishable_colors(size(vorvx{t},2));
% for i = 1:size(vorvx,2)
%         if ~isempty(vorvx{i})
%             k = k+1;
% 
%             plot(vorvx{i}(:,1),vorvx{i}(:,2),'-','Color','b');
%             hold on;
% 
%         end
% end
% plot(pos(:,1),pos(:,2),'o');

%%
% obtain set of voronoi neighbors/vertices
[Abnd,bbnd] = vert2lcon(bnd_pnts,tol);              % obtain series of linear inequalities that defined Voronoi regions
%% create list
for i = 1:size(pos,1)
    k = 0;
    for j = 1:size(vornb{i},2)
        if vornb{i}(1,j) > i
            k = k + 1;
            vornb2{i}(1,k) = vornb{i}(1,j);
        end
    end
end
for m1 =1:size(vornb2,2)
    for j = 1:size(vornb2{m1},2)
        c1 = m1;
        c2 = vornb2{m1}(1,j);
        clear Aag1 bag1 Aag2 bag2 Aagmt1 bagmt1 Aagmt2 bagmt2 pos1 pos2 IDl_tm1 IDl_tm2 Vl_tmp1 Vl_tmp2 voronoi_rg_tmp1 voronoi_rg_tmp2
        % given (c1,c2) where c1< c2
        % remove c2, compute voronoi vertices of c1
        k = 0;
        for i = 1:size(pos,1)
            if i ~= c2
                k = k + 1;
                pos1(k,:) = pos(i,:);
            end
        end
        [~,~,Aag1,bag1] = pVoronoi(pos1,bnd_pnts);
        % obtain intersection with V(c2) from the original Voronoi diagram

        Aagmt1 = [Aag1{c1};A_org{c2};Abnd];
        bagmt1 = [bag1{c1};b_org{c2};bbnd];
        Vl_tmp1= lcon2vert(Aagmt1,bagmt1,[],[],tol);
        if ~isempty(Vl_tmp1)
            % remove Voronoi regions that are not in the interior of the polytope
            if inhull(Vl_tmp1,bnd_pnts,[],1e-15)
%             if inhull(Vl{c1,c2},vorvx{c1},[],0.001) & inhull(Vl{c1,c2},vorvx{c2},[],1.e-13*mean(abs(bnd_pnts(:)))) 
                IDl_tmp1 = convhull(Vl_tmp1);
                voronoi_rg_tmp1 = Vl_tmp1(IDl_tmp1,:);
            else
                voronoi_rg_tmp1 = [];
            end
        else
            voronoi_rg_tmp1 = [];
        end        
        
        % remove c1, compute voronoi vertices of c2
        k = 0;
        for i = 1:size(pos,1)
            if i ~= c1
                k = k +1;
                pos2(k,:) = pos(i,:);
            end
        end
        % obtain intersection with V(c1) from the original Voronoi diagram

        [~,~,Aag2,bag2] = pVoronoi(pos2,bnd_pnts);
        Aagmt2 = [Aag2{c2-1};A_org{c1};Abnd];
        bagmt2 = [bag2{c2-1};b_org{c1};bbnd];
        Vl_tmp2= lcon2vert(Aagmt2,bagmt2,[],[],tol);
        if ~isempty(Vl_tmp2)
            % remove Voronoi regions that are not in the interior of the polytope
            if inhull(Vl_tmp2,bnd_pnts,[],1e-15)
%             if inhull(Vl{c1,c2},vorvx{c1},[],0.001) & inhull(Vl{c1,c2},vorvx{c2},[],1.e-13*mean(abs(bnd_pnts(:)))) 
                IDl_tmp2 = convhull(Vl_tmp2);
                voronoi_rg_tmp2 = Vl_tmp2(IDl_tmp2,:);
            else
                voronoi_rg_tmp2 = [];
            end
        else
            voronoi_rg_tmp2 = [];
        end
        
        voronoi_rg_tmp = [voronoi_rg_tmp1;voronoi_rg_tmp2];
%         [~, idx_i] = unique(num2str(voronoi_rg_tmp,4),'rows'); % precision 4
        if ~isempty(voronoi_rg_tmp ) %&& size(idx_i,1)> size(pos,2)
            voronoi_rg{c1,c2} = voronoi_rg_tmp(convhull(voronoi_rg_tmp),:);
        else
            voronoi_rg{c1,c2} = [];
        end
    end
end


% %%% test
% 
% h0 = figure('position',[0 0 700 700],'Color',[1 1 1]);
% k = 0;
% t = 16;
% for i = 1:size(voronoi_rg,1)*size(voronoi_rg,2)
%     col(i,:)= rand(1,3);
%     col(i,:) = [i/(size(voronoi_rg{t},1)*size(voronoi_rg{t},2)) 1 1];
% end
% col = distinguishable_colors(size(voronoi_rg,1)*size(voronoi_rg,2));
% for i = 1:size(voronoi_rg,1)
%     for j = 1:size(voronoi_rg,2)
%         if ~isempty(voronoi_rg{i,j})
%             k = k+1;
%             if (k == 53)
%                 i
%                 j
%             end
%             plot(voronoi_rg{i,j}(:,1),voronoi_rg{i,j}(:,2),'-','Color','b');
%             hold on;
%             patch(voronoi_rg{i,j}(:,1),voronoi_rg{i,j}(:,2),col(k,:));
%         end
%     end
% end
% bdp = convhull(bnd_pnts);
% plot(bnd_pnts(bdp,1),bnd_pnts(bdp,2),'b-');
% hold on;
% plot(pos(:,1),pos(:,2),'Marker','o','MarkerSize',6,'MarkerFaceColor','r','Color','b','LineStyle','none');hold on;
% plot(p_sav{t}(adv,1),p_sav{t}(adv,2),'Marker','o','MarkerSize',24,'MarkerFaceColor','r','Color','b','LineStyle','none'); hold on;
% axis('equal')
% axis([0 1 0 1]);
% axis('off');
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);  