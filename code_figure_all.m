%% 0.1: load r, two method: ode and perturbation
clear; clc
% % ode
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code\data_use'
load ode_rr_one_dln.mat
num_n=sqrt(size(ode_r_one_dln,1));
ele=zeros(num_n*num_n,3);
ele(:,1)=1:num_n*num_n;
for i=1:num_n
    ele(7*i-6:7*i,2)=i;
    ele(7*i-6:7*i,3)=1:num_n;
end
for i=1:num_n*num_n
    a(i,1)=ismember(0,(ele(i,2)-ele(i,3)));
end
b=find(a==1);
ele_non_dia=ele;ele_non_dia(b,:)=[];
ele_dia=ele(b,:);
a=find(ode_r_one_dln(:,1)~=0);
b=setdiff(a, ele_dia(:,1));
ele_ode_non0=ele(b,:);
ode_r_non0=ode_r_one_dln(b,:);
ode_r_rel_non0=ode_r_non0./(sum(abs(ode_r_non0),2));
clear a b i ode_r_non0
% % data
load data_rr_one_dln.mat
% load data_R_one_dln.mat
% load data_J_one_dln.mat
load data_x.mat
r_one_dln0=r_one_dln;


%% 0.2: calculate related information

% % combian r of E H M at every perturbation
num=2000; % choose N
a_num=zeros(size(r_one_dln0));
for j=1:size(a_num,2)
    for i=1:size(a_num,1)
        a_num(i,j)=size(r_one_dln0{i,j},2);
    end
end
num_min=min([min(a_num),num]);
a_rand=cell(size(r_one_dln0));
for j=1:size(a_rand,2)
    for i=1:size(a_rand,1)
        a_rand{i,j}=randi(a_num(i,j),num_min,1);
    end
end
for j=1:size(r_one_dln0,2)
    for i=1:size(r_one_dln0,1)
        r_one_dln{i,j}=r_one_dln0{i,j}(:,a_rand{i,j});
        % R_one_dln{i,j}=R_one_dln{i,j}(:,a_rand{i,j});
        % J_one_dln{i,j}=J_one_dln{i,j}(:,a_rand{i,j});
        % x{i,j}=x{i,j}(:,:,a_rand{i,j});
    end
end
r_one_dln_2000=r_one_dln;
% save r_one_dln_2000.mat r_one_dln_2000


% % the redefined local response matrix: means and std
% % use the 2000 data r_one_dln_2000 to analysis
% load r_one_dln_2000.mat
% r_one_dln=r_one_dln_2000;
mean_r_one_dln0=cell(size(r_one_dln)); % the redefined local response matrix
std_r_one_dln=cell(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        mean_r_one_dln0{i,j}=mean(r_one_dln{i,j},2);
        std_r_one_dln{i,j}=std(r_one_dln{i,j},[],2);
    end
end
mean_r_one_dln=mean_r_one_dln0;
% % calculate 95% CI
CI=1.96;
CIlr=cell(size(r_one_dln));
CIlrsqrt=cell(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        CIlr{i,j}=[mean_r_one_dln{i,j}-CI*std_r_one_dln{i,j},mean_r_one_dln{i,j}+CI*std_r_one_dln{i,j}];
        id=sign(CIlr{i,j}(:,1)).*sign(CIlr{i,j}(:,2));
        id(id==-1)=0;
        CIlr{i,j}(:,3)=id;
        CIlr{i,j}(:,4)=CIlr{i,j}(:,3).*sign(CIlr{i,j}(:,2));
        CIlrsqrt{i,j}=[mean_r_one_dln{i,j}-CI*std_r_one_dln{i,j}/sqrt(num),...
            mean_r_one_dln{i,j}+CI*std_r_one_dln{i,j}/sqrt(num)];
        id=sign(CIlrsqrt{i,j}(:,1)).*sign(CIlrsqrt{i,j}(:,2));
        id(id==-1)=0;
        CIlrsqrt{i,j}(:,3)=id;
        CIlrsqrt{i,j}(:,4)=CIlrsqrt{i,j}(:,3).*sign(CIlrsqrt{i,j}(:,2));
    end
end
CI_non0=cell(size(r_one_dln));
CIsqrt_non0=cell(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        CI_non0{i,j}=setdiff(find(CIlr{i,j}(:,3)==1),ele_dia(:,1));
        CIsqrt_non0{i,j}=setdiff(find(CIlrsqrt{i,j}(:,3)==1),ele_dia(:,1));
    end
end
% % set rij=0 when 0 in CI
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        id=setdiff(ele_non_dia(:,1),CI_non0{i,j});
        mean_r_one_dln{i,j}(id,:)=0;
    end
end


% % the relative local response matrix
mean_r_rel=cell(size(mean_r_one_dln,1),1); % the relative local response matrix
for i=1:size(mean_r_rel,1)
    a=[mean_r_one_dln{i,1},mean_r_one_dln{i,2},mean_r_one_dln{i,3}];
    mean_r_rel{i,1}=a./sum(abs(a),2);
    mean_r_rel{i,1}(isnan(mean_r_rel{i,1}))=0;
end


% % the non0zero reltive about mean_r and mean_r_rel
mean_r_non0=cell(size(mean_r_one_dln,1),1);
mean_r_rel_non0=cell(size(mean_r_non0));
for i=1:size(mean_r_non0,1)
    id=unique([CI_non0{i,1};CI_non0{i,2};CI_non0{i,3};ele_ode_non0(:,1)]);
    mean_r_non0{i,1}=[mean_r_one_dln{i,1}(id,:),mean_r_one_dln{i,2}(id,:),mean_r_one_dln{i,3}(id,:)];
    mean_r_rel_non0{i,1}=mean_r_rel{i,1}(id,:);
end


% % error about sqrt(mean_r-ode)^2, combian
err_mean_ode_r_dln=zeros(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        err_mean_ode_r_dln(i,j)=sqrt(sum((mean_r_one_dln{i,j}-ode_r_one_dln(:,j))...
            .*(mean_r_one_dln{i,j}-ode_r_one_dln(:,j)),1));
    end
end


% % error and mean std about sqrt(data_r-ode)^2, combian
err_data_ode_r_dln=cell(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        err_data_ode_r_dln{i,j}=sqrt(sum((r_one_dln{i,j}-ode_r_one_dln(:,j))...
            .*(r_one_dln{i,j}-ode_r_one_dln(:,j)),1))';
    end
end
mean_err_data_ode_r_dln=zeros(size(r_one_dln));
std_err_data_ode_r_dln=zeros(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        mean_err_data_ode_r_dln(i,j)=mean(err_data_ode_r_dln{i,j});
        std_err_data_ode_r_dln(i,j)=std(err_data_ode_r_dln{i,j});
    end
end


% % ksdensity pdf
v0=0.025;ro=4;a1=0.01;
r_one_dln0=r_one_dln;
r_one_dln_fre=cell(size(r_one_dln));
r_one_dln_x=cell(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        r_one_dln0{i,j}=round(r_one_dln0{i,j},ro);
        id=find(abs(mean_r_one_dln{i,j})<a1);
        r_one_dln0{i,j}(id,:)=round(r_one_dln0{i,j}(id,:),1);
        for k=1:num_n*num_n
            if all(~(diff(round(r_one_dln0{i,j}(k,:),1))))-1==0
                r_one_dln0{i,j}(k,1)=r_one_dln0{i,j}(k,1)-v0*(mod(k+i,7)+1);
                r_one_dln0{i,j}(k,2)=r_one_dln0{i,j}(k,2)+v0*(mod(k+i,7)+2);
            end
        [r_one_dln_fre{i,j}(k,:),r_one_dln_x{i,j}(k,:)] = ksdensity(r_one_dln0{i,j}(k,:));
        end
    end
end


% % perturbation combian
pt_c=[1; 2; 3; 4]; % perturbation type combian -----------------------------------------------------
pt=length(pt_c);  
p_yla=cell(1,pt+1);
p_yla(1,:)={'','$$\mathrm{I}$$','$$\mathrm{II}$$','$$\mathrm{III}$$','$$\mathrm{IV}$$'};
p_r_x=r_one_dln_x(pt_c,:);
p_r_fre=r_one_dln_fre(pt_c,:);
p_mean_r=mean_r_one_dln(pt_c,:);
p_CI_non0=CI_non0(pt_c,:);
p_mean_r_rel=mean_r_rel(pt_c,:);
p_err_mean_ode_r_dln=err_mean_ode_r_dln(pt_c,:);
max_fre=zeros(num_n*num_n,pt,3);
for k=1:3
    for j=1:pt
        for i=1:num_n*num_n
            max_fre(i,j,k)=max(p_r_fre{j,k}(i,:));
        end
    end
end
clear num a_rand a_num r_one_dln_2000 id i j a k ele_non_dia
clear mean_r_rel mean_r_non0 mean_r_rel_non0 mean_r_one_dln0 std_r_one_dln
clear CI CIlr CIlrsqrt CI_non0 CIsqrt_non0 r_one_dln0 r_one_dln_fre r_one_dln_x v0 ro a1 

%% 0.3: the set befor plot figure
max_r=zeros(size(r_one_dln));
for j=1:size(r_one_dln,2)
    for i=1:size(r_one_dln,1)
        max_r(i,j)=max(mean_r_one_dln{i,j});
    end
end
max_r=max([ max(max(abs(ode_r_one_dln))),max(max(abs(max_r)))]); % the max of r absolute
FS=24; % Fontsize
LW=4; % Linewidth
tit={'E','H','M'};%["E","H","M";"E state","H state","M state"];
xla=cell(size(ele,1),5);
for i=1:size(ele,1)
    xla{i,1}=['$$r_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    xla{i,2}=['$$r^0_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    xla{i,3}=['$$\widetilde{r}^0_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    xla{i,4}=['$$\hat{r}_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    xla{i,5}=['$$\widetilde{r}_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    xla{i,6}=['$$L_{',num2str(ele(i,2)),num2str(ele(i,3)),'}$$'];
    end
path1='E:\huqing\PhD\all_paper\a2_MRA\all_code\figure_all\fig\'; % the path for mat
path2='E:\huqing\PhD\all_paper\a2_MRA\all_code\figure_all\eps\'; % the path for eps
path3='E:\huqing\PhD\all_paper\a2_MRA\all_code\figure_all\jpg\'; % the path for jpg


% acquire the color for plot: 60 sets and six different colors for each set
cd 'E:\huqing\PhD\all_paper\a2_MRA\all_code'
[all_themes, all_colors] = GetColors();

clear i j 


%% figuer 2: pdf of nonzero of 95% CI : E H M 
f=figure;
set(gcf,'Position',[10 60 2700 600])
for ks=1:3
y0=20;
dis=[];
for i=1:pt
    dis=union(dis,p_CI_non0{i,ks});
end
dis=unique(dis);
ytic0=zeros(1,pt+1);ytic0(:,1)=y0;
ytic=zeros(1,pt+1);ytic(:,1)=y0;
id=zeros(size(dis,1),1)+y0/2;
in1=3;
if ks==1
    id(3,1)=id(3,1)+in1;id(4,1)=id(4,1)+in1; % E
    for j=1:pt
        ytic0(1,j+1)=max(max_fre(dis,j,ks));
        ytic(1,j+1)=ytic0(1,j+1)*0.7+ytic(1,j);
    end
elseif ks==2
    id(1,1)=id(1,1)+in1;id(3,1)=id(3,1)+in1;id(6,1)=id(6,1)-in1; % H
    for j=1:pt
        ytic0(1,j+1)=max(max_fre(dis,j,ks));
        ytic(1,j+1)=ytic0(1,j+1)*0.25+ytic(1,j);
    end
elseif ks==3
    id(1,1)=id(1,1)+in1;id(3,1)=id(3,1)+in1; % M
    for j=1:pt
        ytic0(1,j+1)=max(max_fre(dis,j,ks));
        ytic(1,j+1)=ytic0(1,j+1)*0.6+ytic(1,j);
    end
end
if ks==1
di=3;el=48;p_r_x{di,ks}(el,:)=0;p_r_fre{di,ks}(el,:)=0; % E
elseif ks==2
di=2;el=45;p_r_x{di,ks}(el,:)=0;p_r_fre{di,ks}(el,:)=0; % E
di=4;el=45;p_r_x{di,ks}(el,:)=0;p_r_fre{di,ks}(el,:)=0; % E
end
subplot(1,3,ks)
for j=1:pt
    rx=p_r_x{j,ks}(dis,:);% x(i_row,:)=[];
    ry=p_r_fre{j,ks}(dis,:);% y(i_row,:)=[];
    for i=1:length(dis)
        plot(rx(i,:),ry(i,:)+ytic(1,j),'color',all_colors(i,:),'LineWidth',LW-2)
        hold on
        fill([rx(i,:),rx(i,end:-1:1)],[ry(i,:)+ytic(1,j),zeros(1,size(ry,2))+ytic(1,j)],...
            all_colors(i,:),'FaceAlpha',0.55,'EdgeAlpha',0,'EdgeColor','none')
        hold on
        end
    hold on
end
set(gca,'Position',[0.045+0.32*(ks-1),0.11,0.266,0.815],'TickLabelInterpreter','latex');
set(gca, 'Box', 'off', 'FontSize', FS+1, 'LineWidth', LW, 'FontWeight','bold',...
    'XGrid', 'off', 'YGrid', 'on', 'TickDir', 'out', ...
    'XColor', 'k',  'YColor', 'k',...
    'Xlim',[-2.5,1.5],'XTick',[-2.5,-2,-1,0,1,1.5],'XTickLabel',[-2.5,-2,-1,0,1,1.5],...
    'Ylim',[0,ytic(1,end-1)+ytic0(1,end)+10],'YTick',[0,ytic(1,1:end-1)],'YTickLabel',[p_yla(1,:)])%ytic(1,end)+y0
leg=string(zeros(1,length(dis)*2)); leg(:)='';%cell(1,length(dis)*6);
for i=1:length(dis)
    leg(1,2*i)=xla(dis(i),1);%6*i-4
end
hl=legend(leg);
set(hl,'Interpreter','latex','box','off','Location','northwest','Orientation','horizon','NumColumns',7, ...
    'FontSize', FS-2, 'LineWidth', LW-2)
% add top and right line
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % set the axes
set(ax, 'FontSize', FS+1,'LineWidth', LW,'XTick', [],'YTick', []);
end
saveas(f,strcat(path1,['f2','.fig']))
saveas(f,strcat(path3,['f2','.jpg']))
exportgraphics(f,strcat(path2,['f2','.eps']))
clear y0 ytic0 ytic in1 ks di el rx ry leg hl ax i j f dis id


%% figure 3: the specific element of the WT and redefined local response matrix: E, H, M
f=figure;
set(gcf,'Position',[10 60 2700 500])
for ks=1:3
dis=[];
for i=1:pt
    dis=union(dis,p_CI_non0{i,ks});
end
dis=unique(dis);
dis=union(dis,ele_ode_non0(:,1));
dis=unique(dis);
a=zeros(pt+1,size(dis,1)+1);
a(1,1:end-1)=ode_r_one_dln(dis,ks)';
for i=1:pt
    a(i+1,:)=[p_mean_r{i,ks}(dis,1)',p_err_mean_ode_r_dln(i,ks)];
end
subplot(1,3,ks)
SHM=SHeatmap(a(:,1:end-1),'Format','sq');
SHM=SHM.draw();
colormap(slanCM(103)) % the color of numerical value
clim([-round(max_r),round(max_r)]) % the range of the numerical value
SHM.setText('FontSize',FS-10); 
for i=1:size(a,1)
    for j=1:size(a,2)-1
        if a(i,j)==0
            SHM.setTextMN(i,j,'String','0','FontName','Times New Roman','HorizontalAlignment','center')
            SHM.setPatchMN(i,j,'EdgeColor',[1 0.85 0.1],'LineWidth',LW-2) % box
        end
    end
end
set(gca,'Position',[0.044+0.32*(ks-1),0.04,0.27,0.95],'TickLabelInterpreter','latex');
set(gca,'Box', 'on','FontSize', FS-4, 'LineWidth', LW+1, 'FontWeight','bold',...
    'XTickLabel',[xla(dis,4)],'YTickLabel',['WT',p_yla(1,2:end)])%,'YTickLabel',tit
CB=colorbar;
set(CB,'location','eastoutside') 
if ks~= 3
    colorbar('off') % delete colorbar
end
end
saveas(f,strcat(path1,['f3','.fig']))
saveas(f,strcat(path3,['f3','.jpg']))
exportgraphics(f,strcat(path2,['f3','.eps']))
clear dis f ks i a j CB 

%% figure 4: error scatter, the inference error, and the sample error
xli=1:1:3*pt+3;

f=figure;
for i=1:3
    a=zeros(size(err_data_ode_r_dln{1,i},1),pt);
    for j=1:pt
        a(:,j)=err_data_ode_r_dln{pt_c(j),i};
    end
    x=repmat(xli(1+(pt+1)*(i-1):1+(pt+1)*(i-1)+pt-1),num_min,1);
    swarmchart(x,a,5,all_colors((1:pt)+1,:),...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0.2); % scatter of err_data_all
    hold on
    plot(xli(1+(pt+1)*(i-1):1+(pt+1)*(i-1)+pt-1),mean_err_data_ode_r_dln(pt_c,i),...
        'ro','MarkerFaceColor','r','LineWidth', LW-2) % mean of err_data_all
    hold on
    be=errorbar(xli(1+(pt+1)*(i-1):1+(pt+1)*(i-1)+pt-1),mean_err_data_ode_r_dln(pt_c,i),...
        std_err_data_ode_r_dln(pt_c,i)); % std of err_data_all
    set(be, 'LineStyle', 'none', 'Color', 'r', 'LineWidth', LW-2)
    hold on
    plot(xli(1+(pt+1)*(i-1):1+(pt+1)*(i-1)+pt-1),err_mean_ode_r_dln(pt_c,i), ...
        'r+','LineWidth', LW-2) 
    % error of mean_ode
    hold on
end
set(gcf,'Position',[100 100 900 600])
set(gca, 'Box', 'off', 'FontSize', FS+1, 'LineWidth', LW,...   
         'Xlim',[0 3*pt+3],'xtick',[0,0.5*pt+0.5,1.5*pt+1.5,2.5*pt+2.5],'xticklabel',[p_yla(1,1),tit(1,:)], ...
         'XGrid', 'off', 'YGrid', 'on','TickDir', 'out',  ...           
         'XColor', 'k',  'YColor', 'k')  % 0,2.5,7.5,12.5
ylabel('Error','FontSize',FS+2,'FontWeight','bold')  
h=legend(p_yla(1,2:end));
set(h,'Interpreter','latex','Location','north','Orientation','horizon','FontSize', FS-2, 'LineWidth', LW-2)

% add top and right line
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  
set(ax,'LineWidth', LW,'XTick', [],'YTick', []);
saveas(f,strcat(path1,'f4.fig'))
saveas(f,strcat(path3,'f4.jpg'))
exportgraphics(f,strcat(path2,'f4.eps'))

clear xli c x be ax i f af h a


%%  figure 5: the specific element of the redefined local response matrix

v_m=zeros(num_n^2,3);
id=[];
for j=1:3
    v_m(:,j)=p_mean_r{1,j};
    id=union(id,find(v_m(:,j)~=0));
end
id=setdiff(id, ele_dia(:,1));
v_m=v_m(id,:)';
f=figure;
SHM=SHeatmap(round(v_m,4),'Format','sq');
SHM=SHM.draw();
colormap(gca, slanCM(103)) 
clim([-round(max_r),round(max_r)]) 
SHM.setText('FontSize',FS-8); 
for i1=1:size(v_m,1)
    for j=1:size(v_m,2)
        if v_m(i1,j)==0%isnan(v_m(i1,j))%
            SHM.setTextMN(i1,j,'String','0','FontName','Times New Roman','HorizontalAlignment','center')
        end
    end
end
set(gcf,'Position',[50 50 1200 600])
set(gca,'TickLabelInterpreter','latex');
set(gca,'Box', 'on','FontSize', FS-2, 'LineWidth', LW, 'FontWeight','bold', ...
    'XTickLabel',xla(id,4),'YTickLabel',tit(1,:))
CB=colorbar;
CB.Location='northoutside';
set(CB,'YTick',-2:1:2); 
set(CB,'YTickLabel',{'-2','-1','0','1','2'}) 
CB.Position=[0.1298 0.78 0.7754 0.0356];
saveas(f,strcat(path1,'f5.fig'))
saveas(f,strcat(path3,'f5.jpg'))
exportgraphics(f,strcat(path2,'f5.eps'),'Resolution',600)
clear v_m f i1 j CB 


%% figure 7: nonzero of Wt and relative matrix under 10 uniform at E H M state
id1=ele_ode_non0(:,1);
tit1={'WT','$$\mathrm{I}$$'};
f=figure;
set(gcf,'Position',[10 -260 1200 900])
for ks=1:3
a=[ode_r_rel_non0(:,ks),p_mean_r_rel{1,1}(id1,ks)]';
subplot(3,1,ks)
SHM=SHeatmap(round(a,4),'Format','sq');
SHM=SHM.draw();
colormap(slanCM(103)) 
clim([-1,1]) 
SHM.setText('FontSize',FS-8);
i=2;
for j=1:size(a,2)
    if a(i,j)==0
        SHM.setTextMN(i,j,'String','0','FontName','Times New Roman','HorizontalAlignment','center')
    end
    if j==find(id1(:,1)==setdiff(id1,id))
        SHM.setPatchMN(i,j,'EdgeColor',[1 0.85 0.1],'LineWidth',LW-2) 
    end
end
set(gca,'Position',[0.1,0.69-0.31*(ks-1),0.82,0.23],'TickLabelInterpreter','latex');
set(gca,'Box', 'on','FontSize', FS-2, 'LineWidth', LW, 'FontWeight','bold', ...
    'XTickLabel',xla(id1,5),'YTickLabel',tit1(1,:)) % xla(id1,6)
CB=colorbar;
set(CB,'location','northoutside') 
set(CB,'YTick',-1:1/3:1); 
set(CB,'YTickLabel',{'-1','-2/3','-1/3','0','1/3','2/3','1'}) 
if ks~= 1
    colorbar('off') 
end
end
saveas(f,strcat(path1,'f7.fig'))
saveas(f,strcat(path3,'f7.jpg'))
exportgraphics(f,strcat(path2,'f7.eps'),'Resolution',600)
clear f ks CB i j